using Printf: @printf, @sprintf

export read_unk, write_unk

"""
    read_unk(filename)
    read_unk(filename, ::FortranText)
    read_unk(filename, ::FortranBinary)

Read wannier90 `UNK` file for the periodic part of Bloch wavefunctions.

# Return
- `ik`: k-point index, start from 1
- `Ψ`: periodic part of Bloch wavefunctions in real space,
    size = `(n_gx, n_gy, n_gz, n_bands, n_spin)`
"""
function read_unk end

function read_unk(filename::AbstractString, ::FortranText)
    n_spin = occursin("NC", filename) ? 2 : 1

    res = open(filename) do io
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

        return ik, Ψ
    end

    return res
end

function read_unk(filename::AbstractString, ::FortranBinary)
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

function read_unk(filename::AbstractString)
    if isbinary(filename)
        format = FortranBinary()
    else
        format = FortranText()
    end
    ik, Ψ = read_unk(filename, format)

    n_gx, n_gy, n_gz, n_bands, n_spin = size(Ψ)
    @info "Reading unk file" filename ik (n_gx, n_gy, n_gz) n_bands n_spin

    return ik, Ψ
end

"""
    write_unk(filename, ik, Ψ; binary=false)
    write_unk(filename, ik, Ψ, ::FortranText)
    write_unk(filename, ik, Ψ, ::FortranBinary)

Write `UNK` file for the periodic part of Bloch wavefunctions.

# Arguments
- ik: at which kpoint? start from 1
- Ψ: Bloch wavefunctions, `size(Ψ) = (n_gx, n_gy, n_gz, n_bands, n_spin)`

# Keyword arguments
- `binary`: write as Fortran unformatted file
"""
function write_unk end

function write_unk(
    filename::AbstractString, ik::Integer, Ψ::Array{<:Complex,5}, ::FortranText
)
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

function write_unk(
    filename::AbstractString, ik::Integer, Ψ::Array{<:Complex,5}, ::FortranBinary
)
    n_gx, n_gy, n_gz, n_bands, n_spin = size(Ψ)

    # gfortran default integer is 4 bytes
    Tint = Int32

    # not Fortran stream IO, so using `FortranFile`
    io = FortranFile(filename, "w")
    write(io, Tint(n_gx), Tint(n_gy), Tint(n_gz), Tint(ik), Tint(n_bands))

    for ib in 1:n_bands
        for is in 1:n_spin
            write(io, ComplexF64.(Ψ[:, :, :, ib, is]))
        end
    end
    return close(io)
end

function write_unk(filename::AbstractString, ik, Ψ; binary=false)
    n_gx, n_gy, n_gz, n_bands, n_spin = size(Ψ)
    @info "Writing unk file" filename ik (n_gx, n_gy, n_gz) n_bands n_spin

    if binary
        format = FortranBinary()
    else
        format = FortranText()
    end
    return write_unk(filename, ik, Ψ, format)
end
