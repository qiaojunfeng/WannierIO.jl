"""
Simplify the R-vectors of a tight-binding operator, such that the inverse
Fourier transform is a simple sum (and is faster).
"""
abstract type AbstractRvectorReducer end

"""
Wigner-Seitz R-vector reducer, which simplifies the R-vectors by absorbing the
R degeneracies into the operator.

The R-vectors are still the same, the operator is divided by the R-vector degeneracy.

$(TYPEDEF)

# Fields
$(FIELDS)
"""
struct WsRvectorReducer <: AbstractRvectorReducer
    """The new R vectors"""
    Rvectors::Vector{Vec3{Int}}

    """The weights to absorb R degeneracies into operators.
    `weights[iR0]` gives the degeneracy of the original R vector.
    The operator is multiplied by `weights[iR0]` to absorb the degeneracy.
    """
    weights::Vector{Float64}

    function WsRvectorReducer(Rvectors::AbstractVector, weights::AbstractVector)
        length(Rvectors) == length(weights) ||
            throw(ArgumentError("Length of Rvectors must match length of weights"))
        return new(Rvectors, weights)
    end
end

"""
Construct a mapping from the original Wigner-Seitz R-vectors to the new R-vectors, and the weights to absorb the R degeneracies into the operator.

# Arguments
- `Rvectors`: The R-vectors of the original R-space.
- `Rdegens`: The degeneracies of the R-vectors.
"""
function WsRvectorReducer(
    Rvectors::AbstractVector{<:AbstractVector{<:Integer}},
    Rdegens::AbstractVector{<:Integer},
)
    _validate_Rvectors_Rdegens(Rvectors, Rdegens)

    # The R vectors are still the same.
    # Absorb R degeneracies into operator.
    weights = 1 ./ Rdegens

    return WsRvectorReducer(Rvectors, weights)
end

function _validate_Rvectors_Rdegens(Rvectors, Rdegens)
    length(Rvectors) == length(Rdegens) ||
        throw(ArgumentError("Length of Rdegens must match length of Rvectors"))
    all(length.(Rvectors) .== 3) || throw(ArgumentError("All R-vectors must be length 3"))
end

function (reducer::WsRvectorReducer)(operator::AbstractVector{<:AbstractMatrix})
    length(operator) == length(reducer.weights) ||
        throw(ArgumentError("Length of operator must match length of weights"))

    # Absorb R degeneracies into operator.
    operator_new = map(zip(operator, reducer.weights)) do (O, w)
        O .* w
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

    """The index mapping from new R vectors to old R vectors.
    `iR0s = mapping[iR]` gives the mapping from the index `iR` of the new R vector to the indices `iR0s` of the original R vectors.
    """
    mapping::Vector{Vector{Int}}

    """The weights to absorb R and T degeneracies into operators.
    `weights[iR]` gives the weight matrices to absorb the degeneracies of the original R vector and the T vectors for the new R vector.
    The operator is multiplied by `weights[iR]` to absorb the degeneracies.
    """
    weights::Vector{Matrix{Float64}}

    function MdrsRvectorReducer(
        Rvectors::AbstractVector, mapping::AbstractVector, weights::AbstractVector
    )
        length(Rvectors) == length(mapping) ||
            throw(ArgumentError("Length of Rvectors and mapping must match"))
        length(weights) <= length(Rvectors) ||
            throw(ArgumentError("Length of weights must be <= length of Rvectors"))
        return new(Rvectors, mapping, weights)
    end
end

function _validate_Tvectors_Tdegens(Tvectors, Tdegens, nR::Integer)
    nR == length(Tvectors) == length(Tdegens) ||
        throw(ArgumentError("Length of Tvectors and Tdegens must match length of Rvectors"))
    all(length.(Tvectors) .== 3) || throw(ArgumentError("All T-vectors must be length 3"))
    map(Tvectors) do Tvecs
        all(length.(Tvecs) .== 3) || throw(ArgumentError("All T-vectors must be length 3"))
    end
    # Number of Wannier functions
    nw = size(Tvectors[1], 1)
    all(size.(Tvectors) .== (nw, nw)) ||
        throw(ArgumentError("All T-vectors must be of size (nw, nw)"))
    all(size.(Tdegens) .== (nw, nw)) ||
        throw(ArgumentError("All T-degeneracies must be of size (nw, nw)"))
end

function MdrsRvectorReducer(
    Rvectors::AbstractVector{<:AbstractVector{<:Integer}},
    Rdegens::AbstractVector{<:Integer},
    Tvectors::AbstractVector{<:AbstractMatrix{<:AbstractVector{<:Integer}}},
    Tdegens::AbstractVector{<:AbstractMatrix{<:Integer}},
)
    # The length of original R-vectors
    nR0 = length(Rvectors)
    # The number of Wannier functions
    nw = size(Tvectors[1], 1)
    _validate_Rvectors_Rdegens(Rvectors, Rdegens)
    _validate_Tvectors_Tdegens(Tvectors, Tdegens, nR0)

    # Expanded R-vectors
    Rvectors_new = Vector{Vec3{Int}}()
    # TODO maybe keep the same order as the orignal R-vectors in the initial part?
    # Rvectors_new = [vec3(R) for R in Rvectors]
    # Mapping from new R vectors to old R vectors
    mapping = [Int[] for _ in eachindex(Rvectors_new)]
    # Weights to absorb R and T degeneracies into operator.
    # Note the length is the same as the original R vectors, as the weights are
    # multiplied to the original operator.
    weights = [ones(Float64, nw, nw) for _ in 1:nR0]

    # generate expanded R vectors, which contains all the R+T
    for (iR0, R0) in enumerate(Rvectors)
        Tvecs = Tvectors[iR0]
        Nᴿ = Rdegens[iR0]
        # T degeneracy is a matrix
        Mᵀ = Tdegens[iR0]
        weights[iR0] ./= Nᴿ .* Mᵀ
        for n in axes(Mᵀ, 2)
            for m in axes(Mᵀ, 1)
                Nᵀ = Mᵀ[m, n]
                for iT in 1:Nᵀ
                    R0T = vec3(R0 .+ Tvecs[m, n][iT])
                    iR = findfirst(==(R0T), Rvectors_new)
                    if isnothing(iR)
                        push!(Rvectors_new, R0T)
                        push!(mapping, [iR0])
                    else
                        push!(mapping[iR], iR0)
                    end
                end
            end
        end
    end

    return MdrsRvectorReducer(Rvectors_new, mapping, weights)
end

function (reducer::MdrsRvectorReducer)(operator::AbstractVector{<:AbstractMatrix})
    length(operator) == length(reducer.weights) ||
        throw(ArgumentError("Length of operator must match length of weights"))

    op_type = eltype(operator[1])
    op_size = size(operator[1])
    # Simplified operator by absorbing R and T degeneracies
    # With length the same as the new R vectors
    nR = length(reducer.Rvectors)
    return map(1:nR) do iR
        O_new = zeros(op_type, op_size)
        for iR0 in reducer.mapping[iR]
            O = operator[iR0]
            w = reducer.weights[iR0]
            O_new .+= O .* w
        end
        O_new
    end
end
