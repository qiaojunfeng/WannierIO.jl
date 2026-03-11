# XSF format
# Specification from http://www.xcrysden.org/doc/XSF.html

export read_xsf, write_xsf

"""
    $(SIGNATURES)

Read `xsf` file.

# Arguments
- `file`: The name of the input file, or an `IO`.

# Return
- `primvec`: `3 * 3`, Å, each column is a primitive lattice vector
- `convvec`: `3 * 3`, Å, each column is a conventional lattice vector
- `atoms`: `n_atoms` String, atomic labels or numbers
- `atom_positions`: length-`n_atoms` vector, Å, cartesian coordinates
- `origin`: `3`, Å, origin of the grid
- `span_vectors`: `3 * 3`, Å, each column is a spanning vector
- `X`: `nx`, fractional coordinate of grid points along the first spanning vector
- `Y`: `ny`, fractional coordinate of grid points along the second spanning vector
- `Z`: `nz`, fractional coordinate of grid points along the third spanning vector
- `W`: `nx * ny * nz`, volumetric data

!!! note

    Only support reading 1 datagrid in `BLOCK_DATAGRID_3D`.
"""
function read_xsf(io::IO)
    primvec = nothing
    convvec = nothing
    atoms = nothing
    atom_positions = nothing
    origin = nothing
    span_vectors = nothing
    X = Y = Z = nothing
    W = nothing

    block_structure = false
    block_datagrid = false
    while !eof(io)
        line = strip(readline(io))
        if isempty(line) || startswith(line, '#')
            continue
        end

        if occursin("CRYSTAL", line) || occursin("SLAB", line)
            block_structure = true
            block_datagrid = false
            continue
        elseif occursin("BEGIN_BLOCK_DATAGRID_3D", line)
            block_structure = false
            block_datagrid = true
            continue
        end

        if block_structure
            if startswith(line, "PRIMVEC")
                # primitive lattice, each column is a lattice vector
                primvec = zeros(Float64, 3, 3)
                for i in 1:3
                    line = strip(readline(io))
                    primvec[:, i] = parse.(Float64, split(line))
                end
            elseif startswith(line, "CONVVEC")
                # conventional lattice, each column is a lattice vector
                convvec = zeros(Float64, 3, 3)
                for i in 1:3
                    line = strip(readline(io))
                    convvec[:, i] = parse.(Float64, split(line))
                end
            elseif startswith(line, "PRIMCOORD")
                # read atom positions
                line = strip(readline(io))
                n_atom = parse(Int, split(line)[1])
                atoms = Vector{String}(undef, n_atom)
                # each column is a position vector
                atom_positions = zeros(Vec3{Float64}, n_atom)
                for i in 1:n_atom
                    line = strip(readline(io))
                    # might be element label, or atomic number
                    atoms[i] = split(line)[1]
                    atom_positions[i] = Vec3(parse.(Float64, split(line)[2:4])...)
                end
            else
                error("Unexpected line in block_structure: $line")
            end
        elseif block_datagrid
            comment = line  # current line is a comment
            line = strip(readline(io))
            # I only read the 1st data grid, others are ignored
            startswith(line, "BEGIN_DATAGRID_3D") ||
                error("BEGIN_DATAGRID_3D not found in line: $line")
            # identifier = chopprefix(line, "BEGIN_DATAGRID_3D_")
            ngx, ngy, ngz = parse.(Int, split(strip(readline(io))))
            origin = parse.(Float64, split(strip(readline(io))))
            # spanning vectors
            span_vectors = zeros(Float64, 3, 3)
            for i in 1:3
                line = strip(readline(io))
                span_vectors[:, i] = parse.(Float64, split(line))
            end
            # column-major
            W = zeros(Float64, ngx, ngy, ngz)
            idx = 1
            while idx <= ngx * ngy * ngz
                line = split(strip(readline(io)))
                ncol = length(line)
                W[idx:(idx + ncol - 1)] = parse.(Float64, line)
                idx += ncol
            end
            occursin("END_DATAGRID_3D", strip(readline(io))) ||
                error("END_DATAGRID_3D not found in line: $line")
            strip(readline(io)) == "END_BLOCK_DATAGRID_3D" ||
                error("END_BLOCK_DATAGRID_3D not found in line: $line")
            block_datagrid = false
        else
            error("Unexpected block: $line")
        end
    end

    if !isnothing(W)
        n_x, n_y, n_z = size(W)
        X = range(0, 1, n_x)
        Y = range(0, 1, n_y)
        Z = range(0, 1, n_z)
    end

    return (; primvec, convvec, atoms, atom_positions, origin, span_vectors, X, Y, Z, W)
end

function read_xsf(filename::AbstractString)
    return open(filename) do io
        read_xsf(io)
    end
end

"""
    $(SIGNATURES)

Write `xsf` file.

# Arguments
- `file`: The name of the output file, or an `IO`.
- `lattice`: `3 * 3`, Å, each column is a lattice vector
- `atom_positions`: length-`n_atoms` vector, fractional coordinates
- `atom_numbers`: `n_atoms`, atomic numbers
- `origin`: `3`, Å, origin of the grid
- `span_vectors`: `3 * 3`, Å, each column is a spanning vector
- `W`: `nx * ny * nz`, volumetric data
"""
function write_xsf(
    io::IO,
    lattice::AbstractMatrix{T},
    atom_positions::Vector{Vec3{T}},
    atom_numbers::AbstractVector{Int},
    origin::Union{AbstractVector{T},Nothing}=nothing,
    span_vectors::Union{AbstractMatrix{T},Nothing}=nothing,
    W::Union{AbstractArray{T,3},Nothing}=nothing,
) where {T<:Real}
    n_atoms = length(atom_numbers)
    length(atom_positions) == n_atoms || throw(DimensionMismatch("incompatible n_atoms"))
    size(lattice) == (3, 3) || throw(DimensionMismatch("incompatible lattice"))
    isnothing(span_vectors) ||
        size(span_vectors) == (3, 3) ||
        throw(DimensionMismatch("incompatible span_vectors"))

    # header
    @printf(io, "%s\n", default_header())

    @printf(io, "CRYSTAL\n")
    @printf(io, "PRIMVEC\n")
    @printf(io, "%15.10f %15.10f %15.10f\n", lattice[:, 1]...)
    @printf(io, "%15.10f %15.10f %15.10f\n", lattice[:, 2]...)
    @printf(io, "%15.10f %15.10f %15.10f\n", lattice[:, 3]...)
    @printf(io, "CONVVEC\n")
    @printf(io, "%15.10f %15.10f %15.10f\n", lattice[:, 1]...)
    @printf(io, "%15.10f %15.10f %15.10f\n", lattice[:, 2]...)
    @printf(io, "%15.10f %15.10f %15.10f\n", lattice[:, 3]...)
    @printf(io, "PRIMCOORD\n")
    @printf(io, "%d 1\n", n_atoms)
    for i in 1:n_atoms
        pos = lattice * atom_positions[i]
        @printf(io, "%d %15.10f %15.10f %15.10f\n", atom_numbers[i], pos...)
    end

    if !any([isnothing(origin), isnothing(span_vectors), isnothing(W)])
        @printf(io, "\n")
        @printf(io, "BEGIN_BLOCK_DATAGRID_3D\n")
        @printf(io, "3D_field\n")
        @printf(io, "BEGIN_DATAGRID_3D_UNKNOWN\n")

        n_x, n_y, n_z = size(W)
        @printf(io, "%d %d %d\n", n_x, n_y, n_z)
        @printf(io, "%15.10f %15.10f %15.10f\n", origin...)
        @printf(io, "%15.10f %15.10f %15.10f\n", span_vectors[:, 1]...)
        @printf(io, "%15.10f %15.10f %15.10f\n", span_vectors[:, 2]...)
        @printf(io, "%15.10f %15.10f %15.10f\n", span_vectors[:, 3]...)

        # column-major
        ncol = 0
        for k in 1:n_z
            for j in 1:n_y
                for i in 1:n_x
                    @printf(io, " %13.5e", W[i, j, k])
                    ncol += 1
                    if ncol == 6
                        @printf(io, "\n")
                        ncol = 0
                    end
                end
            end
        end
        ncol != 0 && @printf(io, "\n")
        @printf(io, "END_DATAGRID_3D\n")
        @printf(io, "END_BLOCK_DATAGRID_3D\n")
    end

    return nothing
end

function write_xsf(
    filename::AbstractString,
    lattice::AbstractMatrix{T},
    atom_positions::Vector{Vec3{T}},
    atom_numbers::AbstractVector{Int},
    origin::Union{AbstractVector{T},Nothing}=nothing,
    span_vectors::Union{AbstractMatrix{T},Nothing}=nothing,
    W::Union{AbstractArray{T,3},Nothing}=nothing,
) where {T<:Real}
    open(filename, "w") do io
        write_xsf(io, lattice, atom_positions, atom_numbers, origin, span_vectors, W)
    end
    return nothing
end
