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

# Vector{Vector} -> Mat3
Mat3(A::AbstractVector) = Mat3(hcat(A...))

# Mat3 -> Vec3{Vec3}
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
