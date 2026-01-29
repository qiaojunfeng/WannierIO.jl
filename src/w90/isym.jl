export read_isym

"""
A symmetry operation `g = {R|t}` in real space, from the `isym` file.
"""
struct SymOp
    """Comment, usually textual description of the symmetry operation."""
    comment::String

    """Rotation matrix."""
    R::Mat3{Int64}

    """Fractional translation vector."""
    t::Vec3{Float64}

    """Time-reversal flag."""
    time_reversal::Bool

    """SU(2) rotation matrix for spinors."""
    u::SMatrix{2,2,ComplexF64}

    """Index of this symmetry operation in all the symmetry operations."""
    isym::Int64

    """Index of the inverse symmetry operation `g^{-1}`."""
    isym_inv::Int64
end

function Base.show(io::IO, ::MIME"text/plain", s::SymOp)
    print(
        io,
        """SymOp ($(s.comment))
          isym = $(s.isym), isym_inv = $(s.isym_inv)
          R = $(s.R)
          t = $(s.t)
          time_reversal = $(s.time_reversal)
          u = $(s.u)
        """,
    )
end

"""
A representation matrix applied to the Bloch functions for a symmetry operation
from the little group of a k-point.

This is the `d` matrix.

`N` is the number of bands.
"""
struct RepMatBand{N}
    """Index of the IBZ kpoint."""
    ik_ibz::Int64

    """Index of the symmetry operation."""
    isym::Int64

    """Representation matrix acting on the bands."""
    d::SMatrix{N,N,ComplexF64}
end

"""
A representation matrix applied to the Wannier functions for a symmetry
operation in real space.

This is the `D` matrix.

`N` is the number of Wannier functions.
"""
struct RepMatWann{N}
    """Index of the symmetry operation."""
    isym::Int64

    """Representation matrix acting on the Wannier functions."""
    D::SMatrix{N,N,ComplexF64}
end

"""
    $(SIGNATURES)

Read `prefix.isym`.

# Return
A named tuple with the following fields:
- `header::String`: Header line.
- `n_symops::Int64`: Number of symmetry operations.
- `spinors::Bool`: Whether spinors are considered.
- `symops::Vector{SymOp}`: Symmetry operations.
- `nkpts_ibz::Int64`: Number of IBZ kpoints.
- `kpoints_ibz::Vector{Vec3{Float64}}`: IBZ kpoints in fractional coordinates.
- `n_bands::Int64`: Number of bands.
- `n_repmat_band::Int64`: Number of representation matrices for symmetry operations
  in the little groups of all IBZ kpoints.
- `repmat_band::Vector{RepMatBand{n_bands}}`: Representation matrices for symmetry
  operations in the little groups of all IBZ kpoints.
- `n_wann::Int64`: Number of Wannier functions.
- `repmat_wann::Vector{RepMatWann{n_wann}}`: Representation matrices for symmetry operations
  acting on the Wannier functions.
"""
function read_isym(filename::AbstractString)
    return open(filename) do io
        header = readline(io)
        header = strip(header)

        line = split(strip(readline(io)))
        n_symops = parse(Int64, line[1])
        spinors = parse_bool(line[2])

        # Read all symmetry operations
        R = zeros(Int64, 3, 3)
        t = zeros(Float64, 3)
        u = zeros(ComplexF64, 2, 2)
        symops = Vector{SymOp}(undef, n_symops)

        for isym in 1:n_symops
            comment = strip(readline(io))
            for j in 1:3
                line = split(readline(io))
                R[j, :] = parse.(Int64, line)
            end
            line = split(readline(io))
            t .= parse.(Float64, line)
            t_rev = parse_bool(readline(io))
            if spinors
                for j in 1:2, i in 1:2
                    a, b = parse.(Float64, readline(io))
                    u[i, j] = complex(a, b)
                end
            else
                u .= I(2)
            end
            isym_inv = parse(Int64, readline(io))

            symops[isym] = SymOp(comment, R, t, t_rev, u, isym, isym_inv)
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

        # Read little group symmetry operations, the dₘₙ(ĥ, k)
        # Two empty lines
        readline(io)
        readline(io)  # usually is " Representation matrix of G_k"

        # n_repmat_band is the total number of symmetry operations in the
        # little groups of all the IBZ kpoints, ĥ k = k, where k ∈ IBZ
        n_bands, n_repmat_band = parse.(Int64, split(readline(io)))

        repmat_band = Vector{RepMatBand{n_bands}}(undef, n_repmat_band)
        d = zeros(ComplexF64, n_bands, n_bands)

        for irep in 1:n_repmat_band
            ik_ibz, isym, n_elems = parse.(Int64, split(readline(io)))
            # Fill all non-zero elements of the representation matrix
            d .= 0
            for _ in 1:n_elems
                line = split(readline(io))
                m, n = parse.(Int64, line[1:2])
                a, b = parse.(Float64, line[3:4])
                d[m, n] = a + im * b
            end
            repmat_band[irep] = RepMatBand{n_bands}(ik_ibz, isym, d)
        end

        # Read rotation matrix Dₘₙ(ĝ) for Wannier functions
        # Two empty lines
        readline(io)
        readline(io)  # usually is " Rotation matrix of Wannier functions"

        n_wann = parse(Int64, readline(io))

        repmat_wann = Vector{RepMatWann{n_wann}}(undef, n_symops)
        D = zeros(ComplexF64, n_wann, n_wann)

        for _ in 1:n_symops
            isym, n_elems = parse.(Int64, split(readline(io)))
            # Fill all the non-zero elements of the rotation matrix
            D .= 0
            for _ in 1:n_elems
                line = split(readline(io))
                m, n = parse.(Int64, line[1:2])
                a, b = parse.(Float64, line[3:4])
                D[m, n] = a + im * b
            end
            repmat_wann[isym] = RepMatWann{n_wann}(isym, D)
        end

        @info "Read isym file" filename n_symops nkpts_ibz n_bands n_wann
        return (;
            header,
            n_symops,
            spinors,
            symops,
            nkpts_ibz,
            kpoints_ibz,
            n_bands,
            n_repmat_band,
            repmat_band,
            n_wann,
            repmat_wann,
        )
    end
end

"""
    $(SIGNATURES)

Build the index mapping from `ik_ibz` and `isym` to the index in `repmat_band`.
"""
function build_mapping_ik_isym(
    repmat_band::AbstractVector{<:RepMatBand};
    nkpts_ibz::Union{Integer,Nothing}=nothing,
    n_symops::Union{Integer,Nothing}=nothing,
)
    n_repmat = length(repmat_band)
    if isnothing(nkpts_ibz)
        nkpts_ibz = maximum(r.ik_ibz for r in repmat_band)
    end
    if isnothing(n_symops)
        n_symops = maximum(r.isym for r in repmat_band)
    end
    mapping = [Vector{Union{Int64,Nothing}}(nothing, n_symops) for _ in 1:nkpts_ibz]

    for ir in 1:n_repmat
        ik_ibz = repmat_band[ir].ik_ibz
        (0 < ik_ibz <= nkpts_ibz) || error("ik_ibz out of range")
        isym = repmat_band[ir].isym
        (0 < isym <= n_symops) || error("isym out of range")
        mapping[ik_ibz][isym] = ir
    end

    return mapping
end
