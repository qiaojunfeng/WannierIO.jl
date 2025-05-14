using StaticArrays

"""
Length-3 vector type.

For atom positions, kpoints, etc.
"""
const Vec3{T} = SVector{3,T} where {T}

vec3(v::Vec3) = v
vec3(v::AbstractVector) = Vec3(v)

"""
3 x 3 matrix type.

For lattice and reciprocal lattice.
"""
const Mat3{T} = SMatrix{3,3,T,9} where {T}

mat3(A::Mat3) = A
mat3(A::AbstractMatrix) = Mat3(A)

"""
    $(SIGNATURES)

Convert `Vector{Vector}` to `Mat3`. Each vector is a column of the matrix.

!!! note

    This is not defined as a constructor of `Mat3` to avoid type piracy.
"""
mat3(A::AbstractVector) = Mat3(reduce(hcat, A))

"""
    $(SIGNATURES)

Convert `Mat3` to `Vec3{Vec3}`. Each column of the matrix is a vector.

!!! note

    This is not defined as a constructor of `vec3` to avoid type piracy.
"""
vec3(A::Mat3) = Vec3(eachcol(A))

"""
Pair type associating a `Symbol` with a `Vec3`.

Used for win file `atoms_frac` and `kpoint_path`.
"""
const SymbolVec3{T} = Pair{Symbol,Vec3{T}} where {T}

symbolvec3(s, v) = SymbolVec3{eltype(v)}(s, vec3(v))
symbolvec3(s::AbstractString, v) = symbolvec3(Symbol(s), v)
symbolvec3(p::Pair) = symbolvec3(p.first, p.second)
symbolvec3(d::Dict) = symbolvec3(only(d))

abstract type FileFormat end

"""
Fortran formatted IO.
"""
struct FortranText <: FileFormat end

"""
Fortran unformatted IO.
"""
struct FortranBinary <: FileFormat end

"""
Fortran unformatted IO with stream access.

For example, file written using these Fortran code:
```fortran
OPEN(UNIT=11, FILE="ustream.demo", STATUS="NEW", ACCESS="STREAM", FORM="UNFORMATTED")
```
"""
struct FortranBinaryStream <: FileFormat end

"""
Plain text format for Wannier90 `win` and `nnkp` files.

The W90 default `win` or `nnkp` are plain text files but are not
simple arrays of numbers that can be read by `readdlm`, therefore this struct
is used to indicate that the file is plain text but need to be handled
by corresponding functions, e.g., [`read_win`](@ref), [`read_nnkp`](@ref), etc.

This somewhat overlaps with [`FortranText`](@ref), but this one is only
used for small input parameter files e.g. `win` and `nnkp` (in comparison with
the [`Wannier90Toml`](@ref) format), while the [`FortranText`](@ref) (in
comparison with the [`FortranBinary`](@ref) format) is used for large matrices
e.g. `amn`, `mmn`, `eig`, etc.
"""
struct Wannier90Text <: FileFormat end

"""
TOML file format for Wannier90 `win` and `nnkp` files.

Here we introduce a TOML format for `win` and `nnkp`, so that once the `win` or
`nnkp` files are converted into TOML, the TOML files can be loaded by standard
TOML parsers without the headache of writing custom parsers in other Julia packages.

See also [`write_win`](@ref), [`write_nnkp`](@ref), etc.
"""
struct Wannier90Toml <: FileFormat end
