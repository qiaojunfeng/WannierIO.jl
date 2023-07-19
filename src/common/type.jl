using StaticArrays

"""
3 x 3 matrix type.

For lattice and recip_lattice.
"""
const Mat3{T} = SMatrix{3,3,T,9} where {T}

"""
Length-3 vector type.

For atom posistions, kpoints, etc.
"""
const Vec3{T} = SVector{3,T} where {T}

"""
`Vector{Vector}` -> `Mat3`
"""
Mat3(A::AbstractVector) = Mat3(reduce(hcat, A))

"""
`Mat3` -> `Vec3{Vec3}`
"""
Vec3(A::Mat3) = Vec3(eachcol(A))

"""
Pair type associating a `Symbol` with a `Vec3`.

For win file `atoms_frac` and `kpoint_path`.
"""
const SymbolVec3{T} = Pair{Symbol,Vec3{T}} where {T}

SymbolVec3(s, v) = SymbolVec3{eltype(v)}(s, v)
SymbolVec3(s::AbstractString, v) = SymbolVec3(Symbol(s), v)
SymbolVec3(p::Pair) = SymbolVec3(p.first, p.second)
SymbolVec3(d::Dict) = SymbolVec3(only(d))

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
