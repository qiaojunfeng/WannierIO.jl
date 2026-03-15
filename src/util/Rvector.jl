"""
Simplify the R-vectors of a tight-binding operator, such that the inverse
Fourier transform is a simple sum (and is faster).
"""
abstract type AbstractRvectorReducer end

"""
Wigner-Seitz R-vector reducer, which simplifies the R-vectors by absorbing the
R degeneracies into the operator.

The R-vectors are still the same, the operator will be divided by the R-vector
degeneracy.

$(TYPEDEF)

# Fields
$(FIELDS)
"""
struct WsRvectorReducer <: AbstractRvectorReducer
    """The new R vectors (same as the original R vectors)."""
    Rvectors::Vector{Vec3{Int}}

    """The R degeneracies of the (original) R vectors.
    `degens[iR]` gives the degeneracy of the (original) `iR`-th R vector.
    The original `operator[iR]` is divided by `degens[iR]` to absorb the degeneracy.
    """
    degens::Vector{Int}

    function WsRvectorReducer(Rvectors::AbstractVector, degens::AbstractVector)
        _validate_Rvectors_Rdegens(Rvectors, degens)
        return new(Rvectors, degens)
    end
end

function Base.show(io::IO, reducer::WsRvectorReducer)
    print(io, "WsRvectorReducer(n_Rvectors=", length(reducer.Rvectors), ")")
end

function Base.show(io::IO, ::MIME"text/plain", reducer::WsRvectorReducer)
    nR = length(reducer.Rvectors)
    degen_min = nR == 0 ? 0 : minimum(reducer.degens)
    degen_max = nR == 0 ? 0 : maximum(reducer.degens)
    print(
        io,
        """WsRvectorReducer(n_Rvectors = $(nR), degens = Int[$(degen_min), ..., $(degen_max)])""",
    )
end

function _validate_Rvectors_Rdegens(Rvectors, Rdegens)
    length(Rvectors) == length(Rdegens) ||
        throw(ArgumentError("Length of Rdegens must match length of Rvectors"))
    all(length.(Rvectors) .== 3) || throw(ArgumentError("All R-vectors must be length 3"))
end

function (reducer::WsRvectorReducer)(operator::AbstractVector{<:AbstractMatrix})
    length(operator) == length(reducer.degens) ||
        throw(ArgumentError("Length of operator must match length of degens"))

    # Absorb R degeneracies into operator.
    operator_new = map(zip(operator, reducer.degens)) do (O, degen)
        O ./ degen
    end
    return operator_new
end

"""
MDRS R-vector reducer, which simplifies the R-vectors by absorbing the
R degeneracies, T vectors and degeneracies into the operator.

This expands the R-vectors to a new set of R-vectors `R̃ = R + T`,
where `T` is the translation vectors; also divide the operator by the
degeneracies of `R` and `T` vectors.

$(TYPEDEF)

# Fields
$(FIELDS)
"""
struct MdrsRvectorReducer <: AbstractRvectorReducer
    """The new R vectors."""
    Rvectors::Vector{Vec3{Int}}

    """The contribution mapping from new R vectors to original operator entries.
    `mapping[iR]` is a vector of `(iR0, m, n)` tuples, where each tuple means
    one MDRS contribution from `operator[iR0][m, n]` to the new `iR`-th R vector.
    """
    mapping::Vector{Vector{Tuple{Int,Int,Int}}}

    """The combined R and T degeneracies of the original R vectors.
    `degens[iR0]` (`iR0 ∈ iR0s`) gives the degeneracy matrices of the original
    `iR0`-th R vector by combining the original `Rdegens` and the `Tdegens` of
    the T vectors.
    The `operator[iR0]` will be divided by `degens[iR0]` to absorb the degeneracies.
    """
    degens::Vector{Matrix{Int}}

    function MdrsRvectorReducer(
        Rvectors::AbstractVector, mapping::AbstractVector, degens::AbstractVector
    )
        length(Rvectors) == length(mapping) ||
            throw(ArgumentError("Length of Rvectors and mapping must match"))
        nR0 = length(degens)
        nR0 <= length(Rvectors) ||
            throw(ArgumentError("Length of degens must be <= length of Rvectors"))

        # mapping[iR] stores tuples (iR0, m, n), where iR0 indexes
        # original-R arrays (operator/degens).
        if !isempty(mapping)
            has_entry = false
            min_iR0 = typemax(Int)
            max_iR0 = typemin(Int)
            for entries in mapping
                for (iR0, _, _) in entries
                    has_entry = true
                    min_iR0 = min(min_iR0, iR0)
                    max_iR0 = max(max_iR0, iR0)
                end
            end
            if has_entry
                min_iR0 >= 1 || throw(ArgumentError("mapping indices must be >= 1"))
                max_iR0 <= nR0 ||
                    throw(ArgumentError("mapping indices exceed length of degens"))
            end
        end

        return new(Rvectors, mapping, degens)
    end
end

function Base.show(io::IO, reducer::MdrsRvectorReducer)
    print(io, "MdrsRvectorReducer(n_Rvectors=", length(reducer.Rvectors), ")")
end

function Base.show(io::IO, ::MIME"text/plain", reducer::MdrsRvectorReducer)
    nR = length(reducer.Rvectors)
    nR0 = length(reducer.degens)
    nw = nR0 == 0 ? 0 : size(reducer.degens[1], 1)
    print(
        io,
        """MdrsRvectorReducer(n_Rvectors = $(nR), n_Rvectors_original = $(nR0), n_wannier = $(nw))""",
    )
end

function _validate_Tvectors_Tdegens(Tvectors, Tdegens, nR::Integer)
    nR == length(Tvectors) == length(Tdegens) ||
        throw(ArgumentError("Length of Tvectors and Tdegens must match length of Rvectors"))
    nR > 0 || throw(ArgumentError("Length of Rvectors must be > 0"))
    # Number of Wannier functions
    nw = size(Tvectors[1], 1)
    map(Tvectors) do Tmat
        size(Tmat) == (nw, nw) ||
            throw(ArgumentError("All T-vectors must be of size (nw, nw)"))
        for Tvecs in Tmat
            all(length.(Tvecs) .== 3) ||
                throw(ArgumentError("All T-vectors must be length 3"))
        end
    end
    map(Tdegens) do D
        size(D) == (nw, nw) ||
            throw(ArgumentError("All T-degeneracies must be of size (nw, nw)"))
    end
end

function MdrsRvectorReducer(
    Rvectors::AbstractVector{<:AbstractVector{<:Integer}},
    Rdegens::AbstractVector{<:Integer},
    Tvectors::AbstractVector{
        <:AbstractMatrix{<:AbstractVector{<:AbstractVector{<:Integer}}}
    },
    Tdegens::AbstractVector{<:AbstractMatrix{<:Integer}},
)
    # The length of original R-vectors
    nR0 = length(Rvectors)
    # The number of Wannier functions
    nw = size(Tvectors[1], 1)
    _validate_Rvectors_Rdegens(Rvectors, Rdegens)
    _validate_Tvectors_Tdegens(Tvectors, Tdegens, nR0)

    # Expanded R-vectors
    # 1. If empty, the new R+T will be added first.
    # Rvectors_new = Vector{Vec3{Int}}()
    # 2. Otherwise, we can first keep the same order as the orignal R-vectors in
    # the initial part, then the new R+T are added later.
    Rvectors_new = [vec3(R) for R in Rvectors]

    # Mapping from new R vectors to old R-vector matrix entries.
    mapping = [Tuple{Int,Int,Int}[] for _ in eachindex(Rvectors_new)]
    # R and T degeneracies to absorb into operator.
    # Note the length is the same as the original R vectors, as the degeneracies
    # are applied to the original operator.
    degens = [ones(Int, nw, nw) for _ in 1:nR0]

    # generate expanded R vectors, which contains all the R+T
    for (iR0, R0) in enumerate(Rvectors)
        Tvecs = Tvectors[iR0]
        Nᴿ = Rdegens[iR0]
        # T degeneracy is a matrix
        Mᵀ = Tdegens[iR0]
        degens[iR0] .*= Nᴿ .* Mᵀ
        for n in axes(Mᵀ, 2)
            for m in axes(Mᵀ, 1)
                Nᵀ = Mᵀ[m, n]
                for iT in 1:Nᵀ
                    # The new R vector R = R0 + T
                    R = vec3(R0 .+ Tvecs[m, n][iT])
                    iR = findfirst(==(R), Rvectors_new)
                    if isnothing(iR)
                        push!(Rvectors_new, R)
                        push!(mapping, [(iR0, m, n)])
                    else
                        push!(mapping[iR], (iR0, m, n))
                    end
                end
            end
        end
    end

    return MdrsRvectorReducer(Rvectors_new, mapping, degens)
end

function (reducer::MdrsRvectorReducer)(operator::AbstractVector{<:AbstractMatrix})
    length(operator) == length(reducer.degens) ||
        throw(ArgumentError("Length of operator must match length of degens"))

    op_type = eltype(operator[1])
    op_size = size(operator[1])
    # Simplified operator by absorbing R and T degeneracies
    # With length the same as the new R vectors
    nR = length(reducer.Rvectors)
    return map(1:nR) do iR
        O_new = zeros(op_type, op_size)
        for (iR0, m, n) in reducer.mapping[iR]
            O = operator[iR0]
            degen = reducer.degens[iR0]
            O_new[m, n] += O[m, n] / degen[m, n]
        end
        O_new
    end
end
