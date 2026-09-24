"""
    read_u_mat(file)

Read wannier90 `prefix_u.mat` or `prefix_u_dis.mat` file.

# Arguments
- `file`: The name of the input file, or an `IO`.

# Return
- `U`: `Udis` (for disentanglement) or `U` (for maximal localization) matrices,
    array of size `n_bands × n_wann × n_kpts`
- `kpoints`: fractional kpoint coordinates
- `header`: 1st line of the file

!!! warning

    The `wannier90` output `prefix_u_dis.mat` internally sorts the band indices
    according to the disnentanglement window, therefore it can be different from
    the original Bloch states, see the code and comments in [`gauge_matrices_dis`](@ref).
"""
function read_u_mat(io::IO)
    header = String(readstrip(io))
    # for u_dis.mat, nwann <= nbands
    # for u.mat, nbands == nwann
    nkpts, nwann, nbands = parse_vector(readstrip(io), Int)

    kpoints = zeros(Vec3{Float64}, nkpts)
    U = zeros(ComplexF64, nbands, nwann, nkpts)

    for ik in 1:nkpts
        # empty line
        readstrip(io)
        kpoints[ik] = vec3(parse_vector(readstrip(io)))

        for iw in 1:nwann
            for ib in 1:nbands
                vals = parse_vector(readstrip(io))
                U[ib, iw, ik] = vals[1] + im * vals[2]
            end
        end
    end

    return (; U, kpoints, header)
end

function read_u_mat(filename::AbstractString)
    return open(filename) do io
        read_u_mat(io)
    end
end

"""
    write_u_mat(file, U, kpoints; header=default_header(), digits=10)

Write wannier90 `prefix_u.mat` or `prefix_u_dis.mat` file.

# Arguments
- `file`: The name of the output file, or an `IO`.
- `U`: `Udis` (for disentanglement) or `U` (for maximal localization) matrices
- `kpoints`: fractional kpoint coordinates

# Keyword arguments
- `header`: 1st line of the file, optional
- `digits`: decimal digits of every real number; wannier90 writes 10, and
    its free-format reader accepts more. Use 16 to store a `Float64` gauge
    losslessly.

!!! warning

    The `wannier90` output `prefix_u_dis.mat` internally sorts the band indices
    according to the disnentanglement window, therefore it can be different from
    the original Bloch states, see the code and comments in [`gauge_matrices_dis`](@ref).
    This function just writes whatever is inside the input `U` matrix, without
    consider the order of disentanglement window.
"""
function write_u_mat(
        io::IO,
        U::AbstractArray{<:Number, 3},
        kpoints::AbstractVector;
        header::AbstractString = default_header(),
        digits::Integer = 10,
    )
    nbands, nwann, nkpts = size(U)
    nkpts == length(kpoints) || throw(DimensionMismatch("inconsistent number of kpoints"))
    digits >= 1 || throw(ArgumentError("digits must be positive"))
    # wannier90 writes `%15.10f`; the field width grows with the precision so
    # the columns stay aligned
    field = "%$(digits + 5).$(digits)f"
    fmt2 = Printf.Format("  $field  $field\n")
    fmt3 = Printf.Format("  $field  $field  $field\n")

    write(io, header, "\n")
    @printf(io, "%d %d %d\n", nkpts, nwann, nbands)

    for ik in 1:nkpts
        # empty line
        write(io, "\n")
        Printf.format(io, fmt3, kpoints[ik]...)

        for iw in 1:nwann
            for ib in 1:nbands
                u = U[ib, iw, ik]
                Printf.format(io, fmt2, real(u), imag(u))
            end
        end
    end

    return nothing
end

function write_u_mat(
        filename::AbstractString,
        U::AbstractArray{<:Number, 3},
        kpoints::AbstractVector;
        header::AbstractString = default_header(),
        digits::Integer = 10,
    )
    nkpts = size(U, 3)
    nkpts > 0 || throw(ArgumentError("U is empty"))
    nkpts == length(kpoints) || throw(DimensionMismatch("inconsistent number of kpoints"))

    return open(filename, "w") do io
        write_u_mat(io, U, kpoints; header, digits)
    end
end
