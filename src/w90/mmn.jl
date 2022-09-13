using Printf: @printf, @sprintf
using Dates: now

export read_mmn, write_mmn

"""
Read plain text mmn file.
"""
function _read_mmn_fmt(filename::AbstractString)
    io = open(filename)

    header = readline(io)

    line = readline(io)
    n_bands, n_kpts, n_bvecs = map(x -> parse(Int64, x), split(line))

    # overlap matrix
    M = zeros(ComplexF64, n_bands, n_bands, n_bvecs, n_kpts)
    # for each point, list of neighbors, (K) representation
    kpb_k = zeros(Int64, n_bvecs, n_kpts)
    kpb_b = zeros(Int64, 3, n_bvecs, n_kpts)

    while !eof(io)
        for ib in 1:n_bvecs
            line = readline(io)
            arr = split(line)
            k = parse(Int64, arr[1])
            kpb_k[ib, k] = parse(Int64, arr[2])
            kpb_b[:, ib, k] = map(x -> parse(Int64, x), arr[3:5])
            for n in 1:n_bands
                for m in 1:n_bands
                    line = readline(io)
                    arr = split(line)
                    o = parse(Float64, arr[1]) + im * parse(Float64, arr[2])
                    M[m, n, ib, k] = o
                end
            end
        end
    end

    @assert !any(isnan.(M))
    close(io)

    return header, M, kpb_k, kpb_b
end

"""
Read Fortran binary mmn file.
"""
function _read_mmn_bin(filename::AbstractString)
    # I use stream io to write mmn, so I should use plain julia `open`
    io = open(filename)

    header_len = 60
    header = read(io, FString{header_len})

    # gfortran default integer size = 4
    # https://gcc.gnu.org/onlinedocs/gfortran/KIND-Type-Parameters.html
    Tint = Int32
    n_bands = read(io, Tint)
    n_kpts = read(io, Tint)
    n_bvecs = read(io, Tint)

    # overlap matrix
    M = zeros(ComplexF64, n_bands, n_bands, n_bvecs, n_kpts)
    # for each point, list of neighbors, (K) representation
    kpb_k = zeros(Int64, n_bvecs, n_kpts)
    kpb_b = zeros(Int64, 3, n_bvecs, n_kpts)

    while !eof(io)
        for ib in 1:n_bvecs
            k = read(io, Tint)
            kpb_k[ib, k] = read(io, Tint)
            kpb_b[1, ib, k] = read(io, Tint)
            kpb_b[2, ib, k] = read(io, Tint)
            kpb_b[3, ib, k] = read(io, Tint)
            for n in 1:n_bands
                for m in 1:n_bands
                    r = read(io, Float64)
                    i = read(io, Float64)
                    M[m, n, ib, k] = r + im * i
                end
            end
        end
    end

    @assert !any(isnan.(M))
    close(io)

    return header, M, kpb_k, kpb_b
end

"""
    read_mmn(filename::AbstractString)

Read `mmn` file.

Returns a `n_bands * n_bands * n_bvecs * n_kpts` array.
"""
function read_mmn(filename::AbstractString)
    @info "Reading mmn file: $filename"

    if isbinary(filename)
        header, M, kpb_k, kpb_b = _read_mmn_bin(filename)
    else
        header, M, kpb_k, kpb_b = _read_mmn_fmt(filename)
    end

    n_bands, _, n_bvecs, n_kpts = size(M)
    println("  header  = ", header)
    println("  n_bands = ", n_bands)
    println("  n_bvecs = ", n_bvecs)
    println("  n_kpts  = ", n_kpts)
    println()

    return M, kpb_k, kpb_b
end

"""
Write plain text mmn file.
"""
function _write_mmn_fmt(
    filename::AbstractString,
    M::AbstractArray{<:Complex,4},
    kpb_k::AbstractMatrix{Int},
    kpb_b::AbstractArray{Int,3},
    header::AbstractString;
)
    n_bands, _, n_bvecs, n_kpts = size(M)

    open(filename, "w") do io
        header = strip(header)
        write(io, header, "\n")

        @printf(io, "    %d   %d    %d \n", n_bands, n_kpts, n_bvecs)

        for ik in 1:n_kpts
            for ib in 1:n_bvecs
                @printf(io, "%d %d %d %d %d\n", ik, kpb_k[ib, ik], kpb_b[:, ib, ik]...)

                for n in 1:n_bands
                    for m in 1:n_bands
                        o = M[m, n, ib, ik]
                        @printf(io, "  %16.12f  %16.12f\n", real(o), imag(o))
                    end
                end
            end
        end
    end
end

"""
Write binary mmn file.
"""
function _write_mmn_bin(
    filename::AbstractString,
    M::AbstractArray{<:Complex,4},
    kpb_k::AbstractMatrix{Int},
    kpb_b::AbstractArray{Int,3},
    header::AbstractString;
)
    n_bands, _, n_bvecs, n_kpts = size(M)

    # gfortran default integer size = 4
    Tint = Int32
    # I use stream io to write mmn, so I should use plain julia `open`
    open(filename, "w") do io
        header_len = 60
        header = FString(header_len, String(strip(header)))
        write(io, header)

        write(io, Tint(n_bands))
        write(io, Tint(n_kpts))
        write(io, Tint(n_bvecs))

        for ik in 1:n_kpts
            for ib in 1:n_bvecs
                write(io, Tint(ik))
                write(io, Tint(kpb_k[ib, ik]))
                write(io, Tint(kpb_b[1, ib, ik]))
                write(io, Tint(kpb_b[2, ib, ik]))
                write(io, Tint(kpb_b[3, ib, ik]))

                for n in 1:n_bands
                    for m in 1:n_bands
                        o = M[m, n, ib, ik]
                        write(io, Float64(real(o)))
                        write(io, Float64(imag(o)))
                    end
                end
            end
        end
    end
end

"""
    write_mmn(filename, M::Array{ComplexF64,4}, kpb_k, kpb_b, header; binary=false)

Write `mmn` file.

# Arguments
- `filename`: output file name
- `M`: `n_bands * n_bands * n_bvecs * n_kpts` array
- `kpb_k`: `n_bvecs * n_kpts` array
- `kpb_b`: `3 * n_bvecs * n_kpts` array
- `header`: header string

# Keyword arguments
- `binary`: if true write in Fortran binary format
"""
function write_mmn(
    filename::AbstractString,
    M::AbstractArray{<:Complex,4},
    kpb_k::AbstractMatrix{Int},
    kpb_b::AbstractArray{Int,3},
    header::AbstractString;
    binary::Bool=false,
)
    n_bands, _, n_bvecs, n_kpts = size(M)
    n_bands != size(M, 2) && error("M must be n_bands x n_bands x n_bvecs x n_kpts")

    size(kpb_k) != (n_bvecs, n_kpts) && error("kpb_k has wrong size")
    size(kpb_b) != (3, n_bvecs, n_kpts) && error("kpb_b has wrong size")

    if binary
        _write_mmn_bin(filename, M, kpb_k, kpb_b, header)
    else
        _write_mmn_fmt(filename, M, kpb_k, kpb_b, header)
    end

    @info "Written to file: $(filename)"
    println()

    return nothing
end

"""
    write_mmn(filename, M::Array{ComplexF64,4}, kpb_k, kpb_b; binary=false)

Write `mmn` file.

Default header is "Created by WannierIO.jl CURRENT_DATE".

# Arguments
- `filename`: output file name
- `M`: `n_bands * n_bands * n_bvecs * n_kpts` array
- `kpb_k`: `n_bvecs * n_kpts` array
- `kpb_b`: `3 * n_bvecs * n_kpts` array

# Keyword arguments
- `binary`: if true write in Fortran binary format
"""
function write_mmn(
    filename::AbstractString,
    M::AbstractArray{<:Complex,4},
    kpb_k::AbstractMatrix{Int},
    kpb_b::AbstractArray{Int,3};
    binary::Bool=false,
)
    header = @sprintf "Created by WannierIO.jl %s" string(now())
    return write_mmn(filename, M, kpb_k, kpb_b, header; binary=binary)
end
