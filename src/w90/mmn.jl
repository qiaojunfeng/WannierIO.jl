export read_mmn, write_mmn

"""
Container for wannier90 `mmn` data.

$(TYPEDEF)

# Fields

$(FIELDS)
"""
struct Mmn{T <: Real, IT <: Integer}
    "Header line (1st line of the file)"
    header::String

    """Overlap matrices with size `n_bands × n_bands × n_bvecs × n_kpts`."""
    M::Array{Complex{T}, 4}

    """Neighbor kpoint indices with size `n_bvecs × n_kpts`."""
    kpb_k::Matrix{IT}

    """Translation vectors for neighboring kpoints with size `n_bvecs × n_kpts`.

    The translation vector `G` is defined as
    `b = kpoints[kpb_k[ib, ik]] + kpb_G[ib, ik] - kpoints[ik]`,
    where `b` is the `ib`-th bvector of the `ik`-th kpoint.
    """
    kpb_G::Matrix{Vec3{IT}}
end

function Base.show(io::IO, mmn::Mmn)
    n_bands, _, n_bvecs, n_kpts = size(mmn.M)
    return print(io, "Mmn(n_kpts=$(n_kpts), n_bvecs=$(n_bvecs), n_bands=$(n_bands))")
end

function Base.show(io::IO, ::MIME"text/plain", mmn::Mmn)
    n_bands, _, n_bvecs, n_kpts = size(mmn.M)

    return print(
        io,
        """Mmn(
          header: $(mmn.header)
          n_kpts: $(n_kpts)
          n_bvecs: $(n_bvecs)
          n_bands: $(n_bands)
                    M: Array{Complex}($(n_bands)×$(n_bands)×$(n_bvecs)×$(n_kpts))
        )""",
    )
end

"""
    read_mmn(file)
    read_mmn(file, ::FortranText)
    read_mmn(file, ::FortranBinaryStream)

Read wannier90 `mmn` file.

# Arguments
- `file`: The name of the input file, or an `IO`.

# Return
- [`Mmn`](@ref) struct containing the data in the file
"""
function read_mmn end

function read_mmn(io::IO, ::FortranText)
    header = strip(readline(io))

    line = readline(io)
    n_bands, n_kpts, n_bvecs = map(x -> parse(Int, x), split(line))

    # overlap matrix
    M = zeros(ComplexF64, n_bands, n_bands, n_bvecs, n_kpts)
    # for each point, list of neighbors
    kpb_k = zeros(Int, n_bvecs, n_kpts)
    # the translation vector G so that the true bvector
    # b = kpoints[kpb_k] + kpb_G - kpoints[k]
    kpb_G = zeros(Vec3{Int}, n_bvecs, n_kpts)

    while !eof(io)
        for ib in 1:n_bvecs
            line = readline(io)
            arr = split(line)
            ik = parse(Int, arr[1])
            kpb_k[ib, ik] = parse(Int, arr[2])
            kpb_G[ib, ik] = Vec3(parse.(Int, arr[3:5]))
            for n in 1:n_bands
                for m in 1:n_bands
                    line = readline(io)
                    arr = split(line)
                    o = parse(Float64, arr[1]) + im * parse(Float64, arr[2])
                    M[m, n, ib, ik] = o
                end
            end
        end
    end

    all(Mk -> all(Mkb -> !any(isnan.(Mkb)), Mk), M) || error(
        "Some elements in M are NaN, maybe the file is corrupted or not in the correct format",
    )
    return Mmn(String(header), M, kpb_k, kpb_G)
end

function read_mmn(io::IO, ::FortranBinaryStream)
    # I use stream io to write mmn, so I should use plain julia `open`
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
    M = zeros(ComplexF64, n_bands, n_bands, n_bvecs, n_kpts)
    # for each point, list of neighbors, (K) representation
    kpb_k = zeros(Int, n_bvecs, n_kpts)
    kpb_G = zeros(Vec3{Int}, n_bvecs, n_kpts)

    while !eof(io)
        for ib in 1:n_bvecs
            ik = read(io, Tint)
            kpb_k[ib, ik] = read(io, Tint)
            kpb_G[ib, ik] = Vec3{Int}(read(io, Tint), read(io, Tint), read(io, Tint))
            for n in 1:n_bands
                for m in 1:n_bands
                    r = read(io, Float64)
                    i = read(io, Float64)
                    M[m, n, ib, ik] = r + im * i
                end
            end
        end
    end

    all(Mk -> all(Mkb -> !any(isnan.(Mkb)), Mk), M) || error(
        "Some elements in M are NaN, maybe the file is corrupted or not in the correct format",
    )
    return Mmn(String(header), M, kpb_k, kpb_G)
end

function read_mmn(filename::AbstractString, format::AbstractFileFormat)
    return open(filename) do io
        read_mmn(io, format)
    end
end

function read_mmn(file::Union{IO, AbstractString})
    format = detect_fortran_format(file; stream = true)
    mmn = read_mmn(file, format)
    _check_dimensions_M_kpb(mmn.M, mmn.kpb_k, mmn.kpb_G)
    return mmn
end

"""
    write_mmn(file, mmn; binary=false)
    write_mmn(file, mmn, ::FortranText)
    write_mmn(file, mmn, ::FortranBinaryStream)

Write wannier90 `mmn` file.

# Arguments
- `file`: The name of the output file, or an `IO`.
- `mmn`: a [`Mmn`](@ref) struct

# Keyword arguments
- `binary`: if true write in Fortran binary format
"""
function write_mmn end

"""
Check the dimensions between the quantities are consistent.
"""
@inline function _check_dimensions_kpb(kpb_k::AbstractMatrix, kpb_G::AbstractMatrix)
    n_bvecs, n_kpts = size(kpb_k)
    size(kpb_G) == (n_bvecs, n_kpts) || throw(DimensionMismatch("kpb_G has wrong dimensions"))
    all(length.(kpb_G) .== 3) || throw(DimensionMismatch("kpb_G entries are not 3-vectors"))
    return nothing
end

"""
    $(SIGNATURES)
"""
@inline function _check_dimensions_M_kpb(M::AbstractArray, kpb_k::AbstractMatrix, kpb_G::AbstractMatrix)
    _check_dimensions_kpb(kpb_k, kpb_G)

    n_bvecs, n_kpts = size(kpb_k)
    size(M, 4) == n_kpts ||
        throw(DimensionMismatch("M and kpb_k have different n_kpts"))
    size(M, 3) == n_bvecs ||
        throw(DimensionMismatch("M and kpb_k have different n_bvecs"))

    n_bands = size(M, 1)
    n_bands > 0 || throw(ArgumentError("Empty M, n_bands = 0"))
    size(M, 2) == n_bands || throw(DimensionMismatch("M[:, :, ib, ik] are not square matrices"))
    return nothing
end

function write_mmn(io::IO, mmn::Mmn, ::FortranText)
    _check_dimensions_M_kpb(mmn.M, mmn.kpb_k, mmn.kpb_G)
    n_bands, _, n_bvecs, n_kpts = size(mmn.M)

    write(io, strip(mmn.header), "\n")

    @printf(io, "    %d   %d    %d \n", n_bands, n_kpts, n_bvecs)

    for ik in 1:n_kpts
        for ib in 1:n_bvecs
            @printf(io, "%d %d %d %d %d\n", ik, mmn.kpb_k[ib, ik], mmn.kpb_G[ib, ik]...)

            for n in 1:n_bands
                for m in 1:n_bands
                    o = mmn.M[m, n, ib, ik]
                    @printf(io, "  %16.12f  %16.12f\n", real(o), imag(o))
                end
            end
        end
    end
    return nothing
end

function write_mmn(io::IO, mmn::Mmn, ::FortranBinaryStream)
    _check_dimensions_M_kpb(mmn.M, mmn.kpb_k, mmn.kpb_G)
    n_bands, _, n_bvecs, n_kpts = size(mmn.M)

    # gfortran default integer size = 4
    Tint = Int32
    # I use stream io to write mmn, so I should use plain julia `open`
    header_len = 60
    write(io, FString(header_len, String(strip(mmn.header))))

    write(io, Tint(n_bands))
    write(io, Tint(n_kpts))
    write(io, Tint(n_bvecs))

    for ik in 1:n_kpts
        for ib in 1:n_bvecs
            write(io, Tint(ik))
            write(io, Tint(mmn.kpb_k[ib, ik]))
            write(io, Tint(mmn.kpb_G[ib, ik][1]))
            write(io, Tint(mmn.kpb_G[ib, ik][2]))
            write(io, Tint(mmn.kpb_G[ib, ik][3]))

            for n in 1:n_bands
                for m in 1:n_bands
                    o = mmn.M[m, n, ib, ik]
                    write(io, Float64(real(o)))
                    write(io, Float64(imag(o)))
                end
            end
        end
    end
    return nothing
end

function write_mmn(filename::AbstractString, mmn::Mmn, format::AbstractFileFormat)
    open(filename, "w") do io
        write_mmn(io, mmn, format)
    end
    return nothing
end

function write_mmn(filename::AbstractString, mmn::Mmn; binary = false)
    format = fortran_format(; binary, stream = true)
    return write_mmn(filename, mmn, format)
end
