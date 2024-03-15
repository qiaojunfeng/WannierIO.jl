"""
    read_u_mat(filename)

Read wannier90 `prefix_u.mat` or `prefix_u_dis.mat` file.

# Arguments
- `filename`: the input file name

# Return
- `U`: `Udis` (for disentanglement) or `U` (for maximal localization) matrices
- `kpoints`: fractional kpoint coordinates
- `header`: 1st line of the file

!!! warning

    The `wannier90` output `prefix_u_dis.mat` internally sorts the band indices
    according to the disnentanglement window, therefore it can be different from
    the original Bloch states, see the code and comments in [`get_Udis`](@ref).
"""
function read_u_mat(filename::AbstractString)
    return open(filename) do io
        # strip and read line
        srline() = strip(readline(io))

        header = String(srline())
        # for u_dis.mat, nwann <= nbands
        # for u.mat, nbands == nwann
        nkpts, nwann, nbands = parse.(Int, split(srline()))
        @info "Reading u_mat file" filename header nkpts nbands nwann

        kpoints = zeros(Vec3{Float64}, nkpts)
        U = [zeros(ComplexF64, nbands, nwann) for _ in 1:nkpts]

        for ik in 1:nkpts
            # empty line
            srline()
            kpoints[ik] = Vec3(parse.(Float64, split(srline()))...)

            for iw in 1:nwann
                for ib in 1:nbands
                    vals = parse.(Float64, split(srline()))
                    U[ik][ib, iw] = vals[1] + im * vals[2]
                end
            end
        end

        return (; U, kpoints, header)
    end
end

"""
    write_u_mat(filename, U, kpoints; header=default_header())

Write wannier90 `prefix_u.mat` or `prefix_u_dis.mat` file.

# Arguments
- `filename`: the input file name
- `U`: `Udis` (for disentanglement) or `U` (for maximal localization) matrices
- `kpoints`: fractional kpoint coordinates

# Keyword arguments
- `header`: 1st line of the file, optional

!!! warning

    The `wannier90` output `prefix_u_dis.mat` internally sorts the band indices
    according to the disnentanglement window, therefore it can be different from
    the original Bloch states, see the code and comments in [`get_Udis`](@ref).
    This function just writes whatever is inside the input `U` matrix, without
    consider the order of disentanglement window.
"""
function write_u_mat(
    filename::AbstractString,
    U::AbstractVector,
    kpoints::AbstractVector;
    header::AbstractString=default_header(),
)
    nkpts = length(U)
    @assert nkpts > 0 "U is empty"
    @assert nkpts == length(kpoints) "inconsistent number of kpoints"
    nbands, nwann = size(U[1])

    @info "Writing u_mat file" filename header nkpts nbands nwann

    return open(filename, "w") do io
        write(io, header, "\n")
        @printf(io, "%d %d %d\n", nkpts, nwann, nbands)

        for ik in 1:nkpts
            # empty line
            write(io, "\n")
            @printf(io, "  %15.10f  %15.10f  %15.10f\n", kpoints[ik]...)

            for iw in 1:nwann
                for ib in 1:nbands
                    u = U[ik][ib, iw]
                    @printf(io, "  %15.10f  %15.10f\n", real(u), imag(u))
                end
            end
        end
    end
end
