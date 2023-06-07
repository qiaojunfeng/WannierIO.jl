using Printf: @printf, @sprintf

export read_unk, write_unk

"""
Read plain text unk file.
"""
function _read_unk_fmt(filename::AbstractString)
    n_spin = occursin("NC", filename) ? 2 : 1
    io = open(filename)

    line = split(strip(readline(io)))
    n_gx, n_gy, n_gz, ik, n_bands = parse.(Int, line)

    Ψ = zeros(ComplexF64, n_gx, n_gy, n_gz, n_bands, n_spin)

    for ib in 1:n_bands
        for is in 1:n_spin
            for iz in 1:n_gz
                for iy in 1:n_gy
                    for ix in 1:n_gx
                        line = split(strip(readline(io)))
                        v1, v2 = parse.(Float64, line)
                        Ψ[ix, iy, iz, ib, is] = v1 + im * v2
                    end
                end
            end
        end
    end

    close(io)
    return ik, Ψ
end

"""
Read binary unk file.
"""
function _read_unk_bin(filename::AbstractString)
    n_spin = occursin("NC", filename) ? 2 : 1
    # unk files are not in Fortran stream io, so I need to open with FortranFile
    io = FortranFile(filename)

    # gfortran default integer is 4 bytes
    Tint = Int32
    n_gx, n_gy, n_gz, ik, n_bands = read(io, (Tint, 5))

    Ψ = zeros(ComplexF64, n_gx, n_gy, n_gz, n_bands, n_spin)

    for ib in 1:n_bands
        for is in 1:n_spin
            record = Record(io)
            read!(record, view(Ψ, :, :, :, ib, is))
            close(record)
        end
    end

    close(io)
    return ik, Ψ
end

"""
    read_unk(filename::AbstractString)

Read `UNK` file for the periodic part of Bloch wavefunctions.

# Return
- `ik`: k-point index
- `Ψ`: periodic part of Bloch wavefunctions in real space,
    size = `(n_gx, n_gy, n_gz, n_bands, n_spin)`
"""
function read_unk(filename::AbstractString)
    @info "Reading unk file:" filename

    if isbinary(filename)
        ik, Ψ = _read_unk_bin(filename)
    else
        ik, Ψ = _read_unk_fmt(filename)
    end

    n_gx, n_gy, n_gz, n_bands, n_spin = size(Ψ)
    println("  n_gx    = ", n_gx)
    println("  n_gy    = ", n_gy)
    println("  n_gz    = ", n_gz)
    println("  ik      = ", ik)
    println("  n_bands = ", n_bands)
    println("  n_spin  = ", n_spin)
    println()

    # ik: at which kpoint? start from 1
    return ik, Ψ
end

"""
Write plain text unk file.
"""
function _write_unk_fmt(filename::AbstractString, ik::Integer, Ψ::Array{<:Complex,5})
    n_gx, n_gy, n_gz, n_bands, n_spin = size(Ψ)

    open(filename, "w") do io
        @printf(io, " %11d %11d %11d %11d %11d\n", n_gx, n_gy, n_gz, ik, n_bands)

        for ib in 1:n_bands
            for is in 1:n_spin
                for iz in 1:n_gz
                    for iy in 1:n_gy
                        for ix in 1:n_gx
                            v1 = real(Ψ[ix, iy, iz, ib, is])
                            v2 = imag(Ψ[ix, iy, iz, ib, is])
                            @printf(io, " %18.10e %18.10e\n", v1, v2)
                        end
                    end
                end
            end
        end
    end
end

"""
Write binary unk file.
"""
function _write_unk_bin(filename::AbstractString, ik::Integer, Ψ::Array{<:Complex,5})
    n_gx, n_gy, n_gz, n_bands, n_spin = size(Ψ)

    # gfortran default integer is 4 bytes
    Tint = Int32

    # not Fortran stream io, so using `FortranFile`
    io = FortranFile(filename, "w")
    write(io, Tint(n_gx), Tint(n_gy), Tint(n_gz), Tint(ik), Tint(n_bands))

    for ib in 1:n_bands
        for is in 1:n_spin
            write(io, ComplexF64.(Ψ[:, :, :, ib, is]))
        end
    end
    return close(io)
end

"""
    write_unk(filename::AbstractString, ik::Integer, Ψ::Array{Complex,4}; binary=false)

Write `UNK` file for the periodic part of Bloch wavefunctions.

# Arguments
- ik: at which kpoint? start from 1
- Ψ: Bloch wavefunctions, `size(Ψ) = (n_gx, n_gy, n_gz, n_bands, n_spin)`

# Keyword arguments
- `binary`: write as Fortran unformatted file
"""
function write_unk(filename::AbstractString, ik::Integer, Ψ::Array; binary::Bool=false)
    if binary
        _write_unk_bin(filename, ik, Ψ)
    else
        _write_unk_fmt(filename, ik, Ψ)
    end

    n_gx, n_gy, n_gz, n_bands = size(Ψ)
    @info "Written to file: $(filename)"
    println("  n_gx    = ", n_gx)
    println("  n_gy    = ", n_gy)
    println("  n_gz    = ", n_gz)
    println("  ik      = ", ik)
    println("  n_bands = ", n_bands)
    println()

    return nothing
end
