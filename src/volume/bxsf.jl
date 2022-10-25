using Printf: @printf
using Dates: now

# BXSF format
# Specification from http://www.xcrysden.org/doc/XSF.html#__toc__14

export read_bxsf, write_bxsf

"""
    read_bxsf(filename::AbstractString)

Read `bxsf` file.

# Return
- `fermi_energy`: Fermi energy in eV
- `origin`: `3`, Å⁻¹, origin of the grid
- `span_vectors`: `3 * 3`, Å⁻¹, each column is a spanning vector
- `X`: `nx`, fractional coordinate of grid points along the first spanning vector
- `Y`: `ny`, fractional coordinate of grid points along the second spanning vector
- `Z`: `nz`, fractional coordinate of grid points along the third spanning vector
- `E`: `n_bands * nx * ny * nz`, eigenvalues at each grid point
"""
function read_bxsf(filename::AbstractString)
    io = open(filename)

    fermi_energy = nothing
    origin = nothing
    span_vectors = nothing
    X = Y = Z = nothing
    E = nothing

    while !eof(io)
        line = strip(readline(io))
        if isempty(line) || startswith(line, '#')
            continue
        end

        if occursin("BEGIN_INFO", line)
            # skip comments
            line = strip(readline(io))
            while isempty(line) || startswith(line, '#')
                line = strip(readline(io))
            end
            @assert startswith(line, "Fermi Energy:")
            fermi_energy = parse(Float64, split(line, ':')[2])
            line = strip(readline(io))
            @assert line == "END_INFO"
        elseif occursin("BEGIN_BLOCK_BANDGRID_3D", line)
            comment = strip(readline(io))
            line = strip(readline(io))
            # There should be only one data grid
            @assert startswith(line, "BEGIN_BANDGRID_3D")
            # identifier = chopprefix(line, "BEGIN_BANDGRID_3D_")
            n_bands = parse(Int, strip(readline(io)))
            n_x, n_y, n_z = parse.(Int, split(strip(readline(io))))
            origin = parse.(Float64, split(strip(readline(io))))
            # spanning vectors
            span_vectors = zeros(Float64, 3, 3)
            for i in 1:3
                line = strip(readline(io))
                span_vectors[:, i] = parse.(Float64, split(line))
            end
            E = zeros(Float64, n_bands, n_x, n_y, n_z)
            # temp storage for each band, but in row-major
            Eib = similar(E, n_z, n_y, n_x)
            for ib in 1:n_bands
                line = strip(readline(io))
                @assert split(line) == ["BAND:", string(ib)]
                idx = 1
                while idx <= n_x * n_y * n_z
                    line = split(strip(readline(io)))
                    ncol = length(line)
                    Eib[idx:(idx + ncol - 1)] = parse.(Float64, line)
                    idx += ncol
                end
                @assert idx == n_x * n_y * n_z + 1
                # to column-major
                E[ib, :, :, :] = permutedims(Eib, [3, 2, 1])
            end
            @assert occursin("END_BANDGRID_3D", strip(readline(io)))
            @assert strip(readline(io)) == "END_BLOCK_BANDGRID_3D"
        end
    end

    if !isnothing(E)
        _, n_x, n_y, n_z = size(E)
        # the kpoint grid is a general grid, i.e., it includes the last kpoint
        # which is periodic to the first kpoint
        X = range(0, 1, n_x)
        Y = range(0, 1, n_y)
        Z = range(0, 1, n_z)
    end

    return (; fermi_energy, origin, span_vectors, X, Y, Z, E)
end

"""
    write_bxsf(
        filename::AbstractString,
        fermi_energy::T,
        origin::AbstractVector{T},
        span_vectors::AbstractMatrix{T},
        E::AbstractArray{T,4},
    ) where {T<:Real}

Write `bxsf` file.

# Arguments
- `fermi_energy`: Fermi energy in eV
- `origin`: `3`, Å⁻¹, origin of the grid
- `span_vectors`: `3 * 3`, Å⁻¹, each column is a spanning vector
- `E`: `n_bands * nx * ny * nz`, eigenvalues at each grid point
"""
function write_bxsf(
    filename::AbstractString,
    fermi_energy::T,
    origin::AbstractVector{T},
    span_vectors::AbstractMatrix{T},
    E::AbstractArray{T,4},
) where {T<:Real}
    size(origin) == (3,) || error("origin should be a 3-element vector")
    size(span_vectors) == (3, 3) || error("span_vectors should be a 3×3 matrix")

    @info "Writing bxsf file: " filename
    io = open(filename, "w")

    # header
    @printf(io, "BEGIN_INFO\n")
    @printf(io, "  # Created by WannierIO.jl %s\n", string(now()))
    @printf(io, "  Fermi Energy: %21.16f\n", fermi_energy)
    @printf(io, "END_INFO\n\n")

    @printf(io, "BEGIN_BLOCK_BANDGRID_3D\n")
    @printf(io, "from_WannierIO.jl_code\n")
    @printf(io, "BEGIN_BANDGRID_3D_fermi\n")
    n_bands, n_x, n_y, n_z = size(E)
    @printf(io, "%d\n", n_bands)
    @printf(io, "%d %d %d\n", n_x, n_y, n_z)
    @printf(io, "%12.7f %12.7f %12.7f\n", origin...)
    @printf(io, "%12.7f %12.7f %12.7f\n", span_vectors[:, 1]...)
    @printf(io, "%12.7f %12.7f %12.7f\n", span_vectors[:, 2]...)
    @printf(io, "%12.7f %12.7f %12.7f\n", span_vectors[:, 3]...)

    for ib in 1:n_bands
        @printf(io, "BAND: %d\n", ib)
        # row-major
        ncol = 0
        for i in 1:n_x
            for j in 1:n_y
                for k in 1:n_z
                    @printf(io, " %16.8e", E[ib, i, j, k])
                    ncol += 1
                    if ncol == 6
                        @printf(io, "\n")
                        ncol = 0
                    end
                end
            end
        end
        ncol != 0 && @printf(io, "\n")
    end
    @printf(io, "END_BANDGRID_3D\n")
    @printf(io, "END_BLOCK_BANDGRID_3D\n")

    close(io)
    return nothing
end
