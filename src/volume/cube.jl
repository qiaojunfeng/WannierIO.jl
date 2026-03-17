# Cube format
# Specification from http://paulbourke.net/dataformats/cube/

export read_cube, write_cube

"""
Container for `cube` volumetric data.

$(TYPEDEF)

# Fields

$(FIELDS)
"""
struct Cube{T<:Real}
    "Atomic positions in Cartesian coordinates in Å"
    atom_positions::Vector{Vec3{T}}

    "Atomic numbers"
    atom_numbers::Vector{Int}

    "Grid origin in Cartesian coordinates in Å"
    origin::Vec3{T}

    "Voxel vectors, each column is a voxel vector in Å"
    voxel_vectors::Mat3{T}

    "Fractional grid coordinates along first voxel direction"
    X::Vector{T}

    "Fractional grid coordinates along second voxel direction"
    Y::Vector{T}

    "Fractional grid coordinates along third voxel direction"
    Z::Vector{T}

    "Volumetric values, with shape (n_x, n_y, n_z)"
    W::Array{T,3}
end

function Base.show(io::IO, cube::Cube)
    n_atoms = length(cube.atom_positions)
    n_x = length(cube.X)
    n_y = length(cube.Y)
    n_z = length(cube.Z)
    print(io, "Cube(n_atoms=$(n_atoms), grid=$(n_x)×$(n_y)×$(n_z))")
end

function Base.show(io::IO, ::MIME"text/plain", cube::Cube)
    n_atoms = length(cube.atom_positions)
    n_x = length(cube.X)
    n_y = length(cube.Y)
    n_z = length(cube.Z)

    print(
        io,
        """Cube(
          n_atoms: $(n_atoms)
          grid: $(n_x)×$(n_y)×$(n_z)
          origin (Å): $(cube.origin)
          voxel_vectors:
            v₁ = $(cube.voxel_vectors[:, 1])
            v₂ = $(cube.voxel_vectors[:, 2])
            v₃ = $(cube.voxel_vectors[:, 3])
        )""",
    )
end

"""
    $(SIGNATURES)

Read `cube` file.

# Arguments
- `file`: The name of the input file, or an `IO`.

!!! note

    By default, `cube` use Bohr unit, here all returns are in Cartesian coordinates, Å unit.
"""
function read_cube(io::IO)

    # header
    header = readline(io; keep=true)
    header *= readline(io; keep=true)

    line = split(strip(readline(io)))
    n_atoms = parse(Int, line[1])
    origin = vec3(parse_vector(line[2:4])) * Bohr

    # number of voxels in domain
    n_voxels = zeros(Int, 3)
    voxel_vectors = zeros(Float64, 3, 3)
    for i in 1:3
        line = split(strip(readline(io)))
        n_v = parse(Int, line[1])
        n_voxels[i] = n_v
        # bohr unit
        voxel_vectors[:, i] = parse_vector(line[2:4])
    end
    # to Å unit
    voxel_vectors = mat3(voxel_vectors) * Bohr

    atom_positions = Vector{Vec3{Float64}}(undef, n_atoms)
    atom_numbers = zeros(Int, n_atoms)
    for i in 1:n_atoms
        line = split(strip(readline(io)))
        atom_numbers[i] = parse(Int, line[1])
        charge = parse(Float64, line[2])
        # cartesian coordinates, Bohr unit to Å
        atom_positions[i] = vec3(parse_vector(line[3:5])) * Bohr
    end

    n_x, n_y, n_z = n_voxels
    # fractional w.r.t. voxel_vectors
    X = collect(range(0.0, n_x - 1.0, n_x))
    Y = collect(range(0.0, n_y - 1.0, n_y))
    Z = collect(range(0.0, n_z - 1.0, n_z))

    W = zeros(Float64, n_x, n_y, n_z)
    # 6 columns per line
    d, r = divrem(n_z, 6)
    if r > 0
        nline = d + 1
    else
        nline = d
    end
    for ix in 1:n_x
        for iy in 1:n_y
            iz = 1
            for _ in 1:nline
                line = split(strip(readline(io)))
                n = length(line)
                # skip empty line
                if n == 0
                    line = split(strip(readline(io)))
                    n = length(line)
                end
                W[ix, iy, iz:(iz + n - 1)] = parse.(Float64, line)
                iz += n
            end
        end
    end

    return Cube(atom_positions, atom_numbers, origin, voxel_vectors, X, Y, Z, W)
end

function read_cube(filename::AbstractString)
    return open(filename) do io
        read_cube(io)
    end
end

"""
    $(SIGNATURES)

Write `cube` file.

# Arguments
- `file`: The name of the output file, or an `IO`.
- `cube`: a [`Cube`](@ref) struct
"""
function write_cube(io::IO, cube::Cube)
    n_atoms = length(cube.atom_numbers)
    length(cube.atom_positions) == n_atoms ||
        throw(DimensionMismatch("incompatible n_atoms"))
    size(cube.voxel_vectors) == (3, 3) ||
        throw(DimensionMismatch("incompatible voxel_vectors"))
    length(cube.origin) == 3 || throw(DimensionMismatch("origin must be 3-vector"))

    # header
    @printf(io, "%s\n", default_header())
    @printf(io, "outer loop: x, middle loop: y, inner loop: z\n")

    # to Bohr
    origin_bohr = cube.origin ./ Bohr
    @printf(io, "%d %12.6f %12.6f %12.6f\n", n_atoms, origin_bohr...)

    n_xyz = size(cube.W)
    for i in 1:3
        # number of voxels
        n_v = n_xyz[i]
        ax = cube.voxel_vectors[:, i] ./ Bohr
        @printf(io, "%d %12.6f %12.6f %12.6f\n", n_v, ax...)
    end

    for i in 1:n_atoms
        n = cube.atom_numbers[i]
        charge = 1.0
        pos = cube.atom_positions[i] ./ Bohr
        @printf(io, "%d %12.6f %12.6f %12.6f %12.6f\n", n, charge, pos...)
    end

    for ix in 1:n_xyz[1]
        for iy in 1:n_xyz[2]
            for iz in 1:n_xyz[3]
                @printf(io, "%12.6g ", cube.W[ix, iy, iz])

                if (iz % 6 == 0)
                    @printf(io, "\n")
                end
            end
            @printf(io, "\n")
        end
    end

    return nothing
end

function write_cube(filename::AbstractString, cube::Cube)
    open(filename, "w") do io
        write_cube(io, cube)
    end
    return nothing
end
