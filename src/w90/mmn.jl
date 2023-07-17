using Printf: @printf

export read_mmn, write_mmn

"""
    read_mmn(filename)
    read_mmn(filename, ::FortranText)
    read_mmn(filename, ::FortranBinaryStream)

Read wannier90 `mmn` file.

# Return
- `M`: length-`n_kpts` vector, each element is a length-`n_bvecs` vector, then
    each element is a `n_bands * n_bands` matrix
- `kpb_k`: length-`n_kpts` vector, each element is a length-`n_bvecs` vector of
    integers for the indices of the neighboring kpoints
- `kpb_G`: length-`n_kpts` vector, each element is a lenght-`n_bvecs` vector of
    of `Vec3{Int}`, which are the translation vectors
- `header`: 1st line of the file

The translation vector `G` is defined as
`b = kpoints[kpb_k[ik][ib]] + kpb_G[ik][ib] - kpoints[ik]`,
where `b` is the `ib`-th bvector of the `ik`-th kpoint.

The 1st version is a convenience wrapper for the other two.
"""
function read_mmn end

function read_mmn(filename::AbstractString, ::FortranText)
    mmn = open(filename) do io
        header = strip(readline(io))

        line = readline(io)
        n_bands, n_kpts, n_bvecs = map(x -> parse(Int, x), split(line))

        # overlap matrix
        M = [[zeros(ComplexF64, n_bands, n_bands) for _ in 1:n_bvecs] for _ in 1:n_kpts]
        # for each point, list of neighbors
        kpb_k = [zeros(Int, n_bvecs) for _ in 1:n_kpts]
        # the translation vector G so that the true bvector
        # b = kpoints[kpb_k] + kpb_G - kpoints[k]
        kpb_G = [zeros(Vec3{Int}, n_bvecs) for _ in 1:n_kpts]

        while !eof(io)
            for ib in 1:n_bvecs
                line = readline(io)
                arr = split(line)
                ik = parse(Int, arr[1])
                kpb_k[ik][ib] = parse(Int, arr[2])
                kpb_G[ik][ib] = Vec3(parse.(Int, arr[3:5]))
                for n in 1:n_bands
                    for m in 1:n_bands
                        line = readline(io)
                        arr = split(line)
                        o = parse(Float64, arr[1]) + im * parse(Float64, arr[2])
                        M[ik][ib][m, n] = o
                    end
                end
            end
        end

        @assert all(Mk -> all(Mkb -> !any(isnan.(Mkb)), Mk), M)
        return (; M, kpb_k, kpb_G, header)
    end
    return mmn
end

function read_mmn(filename::AbstractString, ::FortranBinaryStream)
    # I use stream io to write mmn, so I should use plain julia `open`
    mmn = open(filename) do io
        header_len = 60
        header = read(io, FString{header_len})
        # from FString to String
        header = strip(String(header))

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
        kpb_G = [zeros(Vec3{Int}, n_bvecs) for _ in 1:n_kpts]

        while !eof(io)
            for ib in 1:n_bvecs
                ik = read(io, Tint)
                kpb_k[ik][ib] = read(io, Tint)
                kpb_G[ik][ib] = Vec3{Int}(read(io, Tint), read(io, Tint), read(io, Tint))
                for n in 1:n_bands
                    for m in 1:n_bands
                        r = read(io, Float64)
                        i = read(io, Float64)
                        M[ik][ib][m, n] = r + im * i
                    end
                end
            end
        end

        @assert all(Mk -> all(Mkb -> !any(isnan.(Mkb)), Mk), M)
        return (; M, kpb_k, kpb_G, header)
    end
    return mmn
end

function read_mmn(filename::AbstractString)
    if isbinary(filename)
        format = FortranBinaryStream()
    else
        format = FortranText()
    end
    M, kpb_k, kpb_G, header = read_mmn(filename, format)

    _check_dimensions_M_kpb(M, kpb_k, kpb_G)
    n_kpts = length(M)
    n_bvecs = length(M[1])
    n_bands = size(M[1][1], 1)
    @info "Reading mmn file" filename header n_kpts n_bvecs n_bands

    # Not returning header since it is printed
    # Note I am returning a Tuple instead of NamedTuple
    return M, kpb_k, kpb_G
end

"""
    write_mmn(filename; M, kpb_k, kpb_G; header=default_header(), binary=false)
    write_mmn(filename; M, kpb_k, kpb_G, ::FortranText; header=default_header(), binary=false)
    write_mmn(filename; M, kpb_k, kpb_G, ::FortranBinaryStream; header=default_header(), binary=false)

Write wannier90 `mmn` file.

# Arguments
- `filename`: output file name
- `M`: length-`n_kpts` vector of `n_bands * n_bands * n_bvecs` arrays
- `kpb_k`: length-`n_kpts` vector of length-`n_bvecs` vector of integers
- `kpb_G`: length-`n_kpts` vector of length-`n_bvecs` vector of `Vec3{Int}` for bvectors

# Keyword arguments
- `header`: header string
- `binary`: if true write in Fortran binary format

The 1st version is a convenience wrapper for the other two.
"""
function write_mmn end

"""
Check the dimensions between the quantities are consistent.
"""
@inline function _check_dimensions_kpb(kpb_k, kpb_G)
    n_kpts = length(kpb_k)
    @assert n_kpts > 0 "Empty kpb_k, n_kpts = 0"
    n_bvecs = length(kpb_k[1])
    @assert n_bvecs > 0 "Empty kpb_k, n_bvecs = 0"

    @assert length(kpb_k) == n_kpts "kpb_k has wrong n_kpts"
    @assert length(kpb_G) == n_kpts "kpb_G has wrong n_kpts"
    @assert all(length.(kpb_k) .== n_bvecs) "kpb_k has different n_bvecs among kpts"
    @assert all(length.(kpb_G) .== n_bvecs) "kpb_G has different n_bvecs among kpts"
    @assert all(all(length.(Gk) .== 3) for Gk in kpb_G) "kpb_G[ib][ib] are not 3-vectors"
end

@inline function _check_dimensions_M_kpb(M, kpb_k, kpb_G)
    _check_dimensions_kpb(kpb_k, kpb_G)

    @assert length(M) == length(kpb_k) "M and kpb_k have different n_kpts"
    n_bvecs = length(kpb_k[1])
    @assert all(length.(M) .== n_bvecs) "M and kpb_k have different n_bvecs"

    n_bands = size(M[1][1], 1)
    @assert n_bands > 0 "Empty M, n_bands = 0"
    @assert all(all(size.(Mk) .== Ref((n_bands, n_bands))) for Mk in M) "M[ik][ib] are not square matrices"
end

function write_mmn(
    filename::AbstractString,
    M::AbstractVector,
    kpb_k::AbstractVector,
    kpb_G::AbstractVector,
    ::FortranText;
    header=default_header(),
)
    _check_dimensions_M_kpb(M, kpb_k, kpb_G)
    n_kpts = length(M)
    n_bvecs = length(M[1])
    n_bands = size(M[1][1], 1)

    open(filename, "w") do io
        header = strip(header)
        write(io, header, "\n")

        @printf(io, "    %d   %d    %d \n", n_bands, n_kpts, n_bvecs)

        for ik in 1:n_kpts
            for ib in 1:n_bvecs
                @printf(io, "%d %d %d %d %d\n", ik, kpb_k[ik][ib], kpb_G[ik][ib]...)

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

function write_mmn(
    filename::AbstractString,
    M::AbstractVector,
    kpb_k::AbstractVector,
    kpb_G::AbstractVector,
    ::FortranBinaryStream;
    header=default_header(),
)
    _check_dimensions_M_kpb(M, kpb_k, kpb_G)
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
            kpb_k_ik = kpb_k[ik]
            kpb_G_ik = kpb_G[ik]
            for ib in 1:n_bvecs
                write(io, Tint(ik))
                write(io, Tint(kpb_k_ik[ib]))
                write(io, Tint(kpb_G_ik[ib][1]))
                write(io, Tint(kpb_G_ik[ib][2]))
                write(io, Tint(kpb_G_ik[ib][3]))

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

function write_mmn(
    filename::AbstractString,
    M::AbstractVector,
    kpb_k::AbstractVector,
    kpb_G::AbstractVector;
    header::AbstractString=default_header(),
    binary::Bool=false,
)
    _check_dimensions_M_kpb(M, kpb_k, kpb_G)
    n_kpts = length(M)
    n_bvecs = length(M[1])
    n_bands = size(M[1][1], 1)

    if binary
        format = FortranBinaryStream()
    else
        format = FortranText()
    end
    write_mmn(filename, M, kpb_k, kpb_G, format; header)

    @info "Writing mmn file" filename header n_kpts n_bvecs n_bands
end
