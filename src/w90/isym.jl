export read_isym, read_isym_raw, standardize

# The `prefix.isym` file is written by Quantum ESPRESSO's pw2wannier90.x
# (`compute_mmn_ibz`, activated by `irr_bz = .true.`). Its symmetry data are
# stored in QE's internal convention, which is *not* the standard Seitz
# convention:
# - the rotation matrix `s` acts on fractional kpoint coordinates as `k' = s k`
# - the pair `(s, ft)` acts on fractional real-space coordinates as
#   `g r = inv(s') (r + ft)`, i.e. `inv(g) r = s' r - ft`
# - the Wannier-function rotation matrices store the representation of the
#   *inverse* operation, `D(inv(g))` (pw2wannier90's `get_rotation_matrix`
#   computes `<g_m|S^-1|g_n>`)
#
# This file therefore provides two layers:
# - a raw layer ([`RawSymOp`](@ref), [`RawIsym`](@ref), [`read_isym_raw`](@ref))
#   that mirrors the file exactly, with fields named after the QE source
#   variables, for debugging and comparison with pw2wannier90
# - a standard layer ([`SymOp`](@ref), [`Isym`](@ref), [`read_isym`](@ref))
#   in the standard (ITA) Seitz convention `g r = W r + v`, produced by
#   [`standardize`](@ref)

"""
A representation matrix ``d(ĥ, k)`` of an element ``ĥ`` of the little group
of an IBZ kpoint, acting on the Bloch states:
```math
ĥ |ψ_{n k}⟩ = \\sum_{n'} |ψ_{n' k}⟩ d_{n' n}(ĥ, k)
```
(the column index is the original state), an `n_bands × n_bands` matrix.
Identical in the raw and standard layers.
"""
struct LittleGroupRep
    """Index of the IBZ kpoint."""
    ik_ibz::Int64

    """Index of the symmetry operation."""
    isym::Int64

    """Representation matrix acting on the Bloch states."""
    d::Matrix{ComplexF64}
end

n_bands(rep::LittleGroupRep) = size(rep.d, 1)
Base.:(==)(a::LittleGroupRep, b::LittleGroupRep) = a.ik_ibz == b.ik_ibz && a.isym == b.isym && a.d == b.d
Base.hash(rep::LittleGroupRep, h::UInt) = hash(rep.d, hash(rep.isym, hash(rep.ik_ibz, hash(:LittleGroupRep, h))))

"""
A representation matrix ``D(ĝ)`` acting on the trial orbitals (and on
symmetry-adapted Wannier functions):
```math
ĝ w^{(0)}_n = \\sum_{n'} D_{n' n}(ĝ) \\, w^{(0)}_{n'}(r - R_{n'}(ĝ))
```
(the column index is the original orbital), an `n_wann × n_wann` matrix.

!!! warning

    In the raw layer ([`RawIsym`](@ref)), entry `isym` stores ``D(ĝ_{isym}^{-1})``,
    the representation of the *inverse* operation, exactly as written by
    pw2wannier90. After [`standardize`](@ref), entry `isym` stores
    ``D(ĝ_{isym})``, its inverse, with the SU(2) sign choice of the `u` field
    of the operation (that of the little-group matrices). An antiunitary
    operation acts as ``D K`` (``K`` complex conjugation).
"""
struct OrbitalRep
    """Index of the symmetry operation."""
    isym::Int64

    """Representation matrix acting on the Wannier functions."""
    D::Matrix{ComplexF64}
end

n_wannier(rep::OrbitalRep) = size(rep.D, 1)
Base.:(==)(a::OrbitalRep, b::OrbitalRep) = a.isym == b.isym && a.D == b.D
Base.hash(rep::OrbitalRep, h::UInt) = hash(rep.D, hash(rep.isym, hash(:OrbitalRep, h)))

"""
A symmetry operation as stored in the `isym` file, in QE's convention.
Field names follow the pw2wannier90 source variables.

Conventions (fractional coordinates):
- kpoints: `k' = s * k` (with `t_rev`: `k' = -s * k`)
- real space: `g r = inv(s') * (r + ft)`
"""
struct RawSymOp
    """Comment, usually textual description of the symmetry operation."""
    comment::String

    """Rotation matrix in QE convention (acts on fractional kpoints)."""
    s::Mat3{Int64}

    """Fractional translation in QE convention."""
    ft::Vec3{Float64}

    """Time-reversal flag."""
    t_rev::Bool

    """SU(2) rotation matrix for spinors."""
    u::SMatrix{2, 2, ComplexF64}

    """Index of this symmetry operation."""
    isym::Int64

    """Index of the inverse symmetry operation."""
    invs::Int64
end

"""
A symmetry operation in the standard (ITA) Seitz convention, all in
fractional coordinates:
- real space: `g r = W * r + v`
- kpoints: `k' = Wk * k` where `Wk = transpose(inv(W))`
  (with `time_reversal`: `k' = -Wk * k`)
"""
struct SymOp
    """Comment, usually textual description of the symmetry operation."""
    comment::String

    """Rotation matrix acting on fractional real-space coordinates (ITA `W`)."""
    W::Mat3{Int64}

    """Fractional translation vector, `g r = W r + v`."""
    v::Vec3{Float64}

    """Cached k-space rotation matrix, `Wk = transpose(inv(W))`, `k' = Wk k`."""
    Wk::Mat3{Int64}

    """Time-reversal flag."""
    time_reversal::Bool

    """SU(2) rotation matrix for spinors."""
    u::SMatrix{2, 2, ComplexF64}

    """Index of this symmetry operation in all the symmetry operations."""
    isym::Int64

    """Index of the inverse symmetry operation `g^{-1}`."""
    isym_inv::Int64
end

function Base.show(io::IO, ::MIME"text/plain", s::SymOp)
    return print(
        io,
        """SymOp ($(s.comment))
          isym = $(s.isym), isym_inv = $(s.isym_inv)
          W = $(s.W)
          v = $(s.v)
          Wk = $(s.Wk)
          time_reversal = $(s.time_reversal)
          u = $(s.u)
        """,
    )
end

"""
Raw container for `prefix.isym` data, mirroring the file exactly.

$(TYPEDEF)

# Fields

$(FIELDS)
"""
struct RawIsym
    "Header line"
    header::String

    "Number of symmetry operations"
    n_symops::Int64

    "Whether spinors are considered"
    spinors::Bool

    "Symmetry operations in QE convention"
    symops::Vector{RawSymOp}

    "Number of IBZ kpoints"
    nkpts_ibz::Int64

    "IBZ kpoints in fractional coordinates"
    kpoints_ibz::Vector{Vec3{Float64}}

    "Number of bands"
    n_bands::Int64

    """Representation matrices `d(ĥ, k)` for the little groups of all IBZ
    kpoints (file section `Representation matrix of G_k`)"""
    repmat_band::Vector{LittleGroupRep}

    "Number of Wannier functions"
    n_wann::Int64

    """Representation matrices for the Wannier functions, storing `D(inv(g))`
    at index `isym` (file section `Rotation matrix of Wannier functions`)"""
    repmat_wann::Vector{OrbitalRep}
end

"""
Container for `prefix.isym` data in the standard Seitz convention.
Produced by [`standardize`](@ref) (or directly by [`read_isym`](@ref)).

$(TYPEDEF)

# Fields

$(FIELDS)
"""
struct Isym
    "Header line"
    header::String

    "Number of symmetry operations"
    n_symops::Int64

    "Whether spinors are considered"
    spinors::Bool

    "Symmetry operations in standard Seitz convention"
    symops::Vector{SymOp}

    "Number of IBZ kpoints"
    nkpts_ibz::Int64

    "IBZ kpoints in fractional coordinates"
    kpoints_ibz::Vector{Vec3{Float64}}

    "Number of bands"
    n_bands::Int64

    """Little-group representation matrices `d(ĥ, k)`, sparse in
    `(ik_ibz, isym)`: only elements of the little group of each IBZ kpoint
    are present"""
    littlegroup_reps::Vector{LittleGroupRep}

    "Number of Wannier functions"
    n_wann::Int64

    """Orbital representation matrices, `orbital_reps[isym]` stores
    `D(g_isym)`, dense in `isym`"""
    orbital_reps::Vector{OrbitalRep}
end

n_bands(x::Union{RawIsym, Isym}) = x.n_bands
n_wannier(x::Union{RawIsym, Isym}) = x.n_wann

"""
    $(SIGNATURES)

Read `prefix.isym` without any convention conversion.

The returned [`RawIsym`](@ref) mirrors the file (QE conventions); use it for
debugging and direct comparison with pw2wannier90. For everything else use
[`read_isym`](@ref), which returns the standard Seitz convention.
"""
function read_isym_raw(io::IO)
    header = readline(io)
    header = strip(header)

    line = split(strip(readline(io)))
    n_symops = parse(Int64, line[1])
    spinors = parse_bool(line[2])

    # Read all symmetry operations
    s = zeros(Int64, 3, 3)
    ft = zeros(Float64, 3)
    u = zeros(ComplexF64, 2, 2)
    symops = Vector{RawSymOp}(undef, n_symops)

    for isym in 1:n_symops
        comment = strip(readline(io))
        for j in 1:3
            line = split(readline(io))
            s[j, :] = parse.(Int64, line)
        end
        line = split(readline(io))
        ft .= parse.(Float64, line)
        t_rev = parse_bool(readline(io))
        if spinors
            for i in 1:2, j in 1:2
                a, b = parse.(Float64, split(readline(io)))
                u[i, j] = complex(a, b)
            end
        else
            u .= I(2)
        end
        invs = parse(Int64, readline(io))

        symops[isym] = RawSymOp(comment, s, ft, t_rev, u, isym, invs)
    end

    # Read IBZ kpoints
    # Two empty lines
    readline(io)
    readline(io)  # usually is " K points"

    nkpts_ibz = parse(Int64, strip(readline(io)))
    # IBZ kpoints in fractional coordinates
    kpoints_ibz = Vector{Vec3{Float64}}(undef, nkpts_ibz)

    for ik in 1:nkpts_ibz
        kpoints_ibz[ik] = parse.(Float64, split(readline(io)))
    end

    # Read little group symmetry operations, the dₘₙ(ĥ, k)
    # Two empty lines
    readline(io)
    readline(io)  # usually is " Representation matrix of G_k"

    # n_repmat_band is the total number of symmetry operations in the
    # little groups of all the IBZ kpoints, ĥ k = k, where k ∈ IBZ
    n_bands, n_repmat_band = parse.(Int64, split(readline(io)))

    repmat_band = Vector{LittleGroupRep}(undef, n_repmat_band)

    for irep in 1:n_repmat_band
        ik_ibz, isym, n_elems = parse.(Int64, split(readline(io)))
        # Fill all non-zero elements of the representation matrix
        d = zeros(ComplexF64, n_bands, n_bands)
        for _ in 1:n_elems
            line = split(readline(io))
            m, n = parse.(Int64, line[1:2])
            a, b = parse.(Float64, line[3:4])
            d[m, n] = a + im * b
        end
        repmat_band[irep] = LittleGroupRep(ik_ibz, isym, d)
    end

    # Read rotation matrix Dₘₙ(ĝ⁻¹) for Wannier functions
    # Two empty lines
    readline(io)
    readline(io)  # usually is " Rotation matrix of Wannier functions"

    n_wann = parse(Int64, readline(io))

    repmat_wann = Vector{OrbitalRep}(undef, n_symops)

    for i in 1:n_symops
        isym, n_elems = parse.(Int64, split(readline(io)))
        # Fill all the non-zero elements of the rotation matrix
        D = zeros(ComplexF64, n_wann, n_wann)
        for _ in 1:n_elems
            line = split(readline(io))
            m, n = parse.(Int64, line[1:2])
            a, b = parse.(Float64, line[3:4])
            D[m, n] = a + im * b
        end
        repmat_wann[i] = OrbitalRep(isym, D)
    end

    return RawIsym(
        String(header),
        n_symops,
        spinors,
        symops,
        nkpts_ibz,
        kpoints_ibz,
        n_bands,
        repmat_band,
        n_wann,
        repmat_wann,
    )
end

function read_isym_raw(filename::AbstractString)
    return open(filename) do io
        read_isym_raw(io)
    end
end

"""
    $(SIGNATURES)

Convert a [`RawSymOp`](@ref) (QE convention) to a [`SymOp`](@ref) (standard
Seitz convention).

The conversion is
```
W  = transpose(inv(s))    # real-space rotation, g r = W r + v
v  = W * ft               # translation
Wk = s                    # k-space rotation, k' = Wk k
```
"""
function standardize(op::RawSymOp)
    Winv_float = transpose(inv(Matrix{Float64}(op.s)))
    W = Mat3{Int64}(round.(Int64, Winv_float))
    isapprox(W, Winv_float; atol = 1.0e-8) ||
        error("transpose(inv(s)) is not an integer matrix for isym = $(op.isym)")
    v = Vec3{Float64}(W * op.ft)
    return SymOp(op.comment, W, v, op.s, op.t_rev, op.u, op.isym, op.invs)
end

"""
    $(SIGNATURES)

Convert a [`RawIsym`](@ref) (file/QE conventions) to an [`Isym`](@ref) in the
standard Seitz convention:
- symmetry operations are converted by [`standardize(::RawSymOp)`](@ref)
- `orbital_reps[isym]` stores `D(g_isym)`, the inverse of the raw entry
  `repmat_wann[isym]`, which is the matrix of the exact inverse operation
  `g_isym^{-1}` (`D^†` for a unitary operation; for an antiunitary one, which
  acts as `D K` with `K` complex conjugation, the inverse acts as `D^T K`,
  so `D` is the transpose of the raw entry)
- `littlegroup_reps` are copied unchanged (the file's `d` matrices are
  already the standard `d(ĥ, k) = ⟨ψ_m|ĥ ψ_n⟩`)

For spinors all three objects then share one SU(2) sign choice (`u(g)` or
`-u(g)` for each operation), that of the `u` field of the operations:
pw2wannier90 rotates the wave functions with `u(g)` and builds the raw
orbital entry of `g^{-1}` with `u(g)^†`. Re-indexing the raw entry of the
stored inverse element `invs(isym)` instead would flip the sign for some
operations (the twofold rotations and mirrors), because the SU(2) matrix of
the stored inverse is `±u(g)^{-1}`.
"""
function standardize(raw::RawIsym)
    symops = standardize.(raw.symops)

    isym2entry = Dict(rep.isym => rep for rep in raw.repmat_wann)
    orbital_reps = map(1:raw.n_symops) do isym
        rep = isym2entry[isym]
        D = raw.symops[isym].t_rev ? transpose(rep.D) : adjoint(rep.D)
        OrbitalRep(isym, Matrix(D))
    end

    return Isym(
        raw.header,
        raw.n_symops,
        raw.spinors,
        symops,
        raw.nkpts_ibz,
        raw.kpoints_ibz,
        raw.n_bands,
        raw.repmat_band,
        raw.n_wann,
        orbital_reps,
    )
end

"""
    $(SIGNATURES)

Read `prefix.isym` and convert to the standard Seitz convention.

Equivalent to `standardize(read_isym_raw(filename))`; see
[`standardize`](@ref) for the conversion and [`Isym`](@ref) for the
conventions of the returned struct.
"""
read_isym(io_or_filename) = standardize(read_isym_raw(io_or_filename))

"""
    $(SIGNATURES)

Build the index mapping from `ik_ibz` and `isym` to the index in
`littlegroup_reps`.
"""
function build_mapping_ik_isym(
        littlegroup_reps::AbstractVector{<:LittleGroupRep};
        nkpts_ibz::Union{Integer, Nothing} = nothing,
        n_symops::Union{Integer, Nothing} = nothing,
    )
    n_reps = length(littlegroup_reps)
    if isnothing(nkpts_ibz)
        nkpts_ibz = maximum(r.ik_ibz for r in littlegroup_reps)
    end
    if isnothing(n_symops)
        n_symops = maximum(r.isym for r in littlegroup_reps)
    end
    mapping = [Vector{Union{Int64, Nothing}}(nothing, n_symops) for _ in 1:nkpts_ibz]

    for ir in 1:n_reps
        ik_ibz = littlegroup_reps[ir].ik_ibz
        (0 < ik_ibz <= nkpts_ibz) || throw(ArgumentError("ik_ibz out of range"))
        isym = littlegroup_reps[ir].isym
        (0 < isym <= n_symops) || throw(ArgumentError("isym out of range"))
        mapping[ik_ibz][isym] = ir
    end

    return mapping
end
