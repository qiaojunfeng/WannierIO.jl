using Printf: @printf

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
    M = [[zeros(ComplexF64, n_bands, n_bands) for _ in 1:n_bvecs] for _ in 1:n_kpts]
    # for each point, list of neighbors, (K) representation
    kpb_k = [zeros(Int, n_bvecs) for _ in 1:n_kpts]
    kpb_b = [zeros(Vec3{Int}, n_bvecs) for _ in 1:n_kpts]

    while !eof(io)
        for ib in 1:n_bvecs
            line = readline(io)
            arr = split(line)
            k = parse(Int64, arr[1])
            kpb_k[k][ib] = parse(Int64, arr[2])
            kpb_b[k][ib] = Vec3(parse.(Int, arr[3:5]))
            for n in 1:n_bands
                for m in 1:n_bands
                    line = readline(io)
                    arr = split(line)
                    o = parse(Float64, arr[1]) + im * parse(Float64, arr[2])
                    M[k][ib][m, n] = o
                end
            end
        end
    end

    @assert all(Mk -> all(Mkb -> !any(isnan.(Mkb)), Mk), M)
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
    M = [[zeros(ComplexF64, n_bands, n_bands) for _ in 1:n_bvecs] for _ in 1:n_kpts]
    # for each point, list of neighbors, (K) representation
    kpb_k = [zeros(Int, n_bvecs) for _ in 1:n_kpts]
    kpb_b = [zeros(Vec3{Int}, n_bvecs) for _ in 1:n_kpts]

    while !eof(io)
        for ib in 1:n_bvecs
            k = read(io, Tint)
            kpb_k[k][ib] = read(io, Tint)
            kpb_b[k][ib] = Vec3{Int}(read(io, Tint), read(io, Tint), read(io, Tint))
            for n in 1:n_bands
                for m in 1:n_bands
                    r = read(io, Float64)
                    i = read(io, Float64)
                    M[k][ib][m, n] = r + im * i
                end
            end
        end
    end

    @assert all(Mk -> all(Mkb -> !any(isnan.(Mkb)), Mk), M)
    close(io)

    return header, M, kpb_k, kpb_b
end

"""
    read_mmn(filename::AbstractString)

Read `mmn` file.

# Return
- `M`: length-`n_kpts` vector of `n_bands * n_bands * n_bvecs` arrays
- `kpb_k`: length-`n_kpts` vector of length-`n_bvecs` vector of integers
- `kpb_b`: length-`n_kpts` vector of length-`n_bvecs` vector of `Vec3{Int}` for bvectors
"""
function read_mmn(filename::AbstractString)
    @info "Reading mmn file: $filename"

    if isbinary(filename)
        header, M, kpb_k, kpb_b = _read_mmn_bin(filename)
    else
        header, M, kpb_k, kpb_b = _read_mmn_fmt(filename)
    end

    n_kpts = length(M)
    n_bvecs = length(M[1])
    n_bands = size(M[1][1], 1)
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
    M::AbstractVector,
    kpb_k::AbstractVector,
    kpb_b::AbstractVector,
    header::AbstractString;
)
    n_kpts = length(M)
    n_bvecs = length(M[1])
    n_bands = size(M[1][1], 1)

    open(filename, "w") do io
        header = strip(header)
        write(io, header, "\n")

        @printf(io, "    %d   %d    %d \n", n_bands, n_kpts, n_bvecs)

        for ik in 1:n_kpts
            for ib in 1:n_bvecs
                @printf(io, "%d %d %d %d %d\n", ik, kpb_k[ik][ib], kpb_b[ik][ib]...)

                for n in 1:n_bands
                    for m in 1:n_bands
                        o = M[ik][ib][m, n]
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
    M::AbstractVector,
    kpb_k::AbstractVector,
    kpb_b::AbstractVector,
    header::AbstractString;
)
    n_kpts = length(M)
    n_bvecs = length(M[1])
    n_bands = size(M[1][1], 1)

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
            kpbk = kpb_k[ik]
            kpbb = kpb_b[ik]
            for ib in 1:n_bvecs
                write(io, Tint(ik))
                write(io, Tint(kpbk[ib]))
                write(io, Tint(kpbb[ib][1]))
                write(io, Tint(kpbb[ib][2]))
                write(io, Tint(kpbb[ib][3]))

                for n in 1:n_bands
                    for m in 1:n_bands
                        o = M[ik][ib][m, n]
                        write(io, Float64(real(o)))
                        write(io, Float64(imag(o)))
                    end
                end
            end
        end
    end
end

"""
    write_mmn(filename, M, kpb_k, kpb_b, header; binary=false)

Write `mmn` file.

# Arguments
- `filename`: output file name
- `M`: length-`n_kpts` vector of `n_bands * n_bands * n_bvecs` arrays
- `kpb_k`: length-`n_kpts` vector of length-`n_bvecs` vector of integers
- `kpb_b`: length-`n_kpts` vector of length-`n_bvecs` vector of `Vec3{Int}` for bvectors
- `header`: header string

# Keyword arguments
- `binary`: if true write in Fortran binary format
"""
function write_mmn(
    filename::AbstractString,
    M::AbstractVector,
    kpb_k::AbstractVector,
    kpb_b::AbstractVector,
    header::AbstractString=default_header();
    binary::Bool=false,
)
    n_kpts = length(M)
    n_bvecs = length(M[1])
    n_bands = size(M[1][1], 1)
    n_bands != size(M[1][1], 2) &&
        error("M at each kpoint and bvector must be n_bands x n_bands")

    (length(kpb_k), length(kpb_k[1])) != (n_kpts, n_bvecs) && error("kpb_k has wrong size")
    (length(kpb_b), length(kpb_b[1])) != (n_kpts, n_bvecs) && error("kpb_b has wrong size")

    if binary
        _write_mmn_bin(filename, M, kpb_k, kpb_b, header)
    else
        _write_mmn_fmt(filename, M, kpb_k, kpb_b, header)
    end

    @info "Written to file: $(filename)"
    println()

    return nothing
end
