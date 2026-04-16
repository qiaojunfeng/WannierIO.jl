export RvectorReducer

"""
Simplify the R-vectors of a tight-binding operator, such that the inverse
Fourier transform is a simple sum (and is faster).
"""
abstract type AbstractRvectorReducer end

"""
Reduce the R-vectors of a tight-binding operator by absorbing the R degeneracies
(and T vectors & degeneracies, if provided) into the operator.

This is a general interface for R-vector reduction, which can be implemented for
different input argument types.
"""
function RvectorReducer end

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

function _validate_Rvectors_Rdegens(Rvectors::AbstractVector, Rdegens::AbstractVector)
    length(Rvectors) == length(Rdegens) ||
        throw(ArgumentError("Length of Rdegens must match length of Rvectors"))
    return all(length.(Rvectors) .== 3) || throw(ArgumentError("All R-vectors must be length 3"))
end

function Base.show(io::IO, reducer::WsRvectorReducer)
    return print(io, "WsRvectorReducer(n_Rvecs=", length(reducer.Rvectors), ")")
end

function Base.:(==)(a::WsRvectorReducer, b::WsRvectorReducer)
    return a.Rvectors == b.Rvectors && a.degens == b.degens
end

function Base.isequal(a::WsRvectorReducer, b::WsRvectorReducer)
    return isequal(a.Rvectors, b.Rvectors) && isequal(a.degens, b.degens)
end

function Base.show(io::IO, ::MIME"text/plain", reducer::WsRvectorReducer)
    nR = length(reducer.Rvectors)
    degen_min = nR == 0 ? 0 : minimum(reducer.degens)
    degen_max = nR == 0 ? 0 : maximum(reducer.degens)
    return print(
        io,
        """WsRvectorReducer(n_Rvecs = $(nR), degens = Int[$(degen_min), ..., $(degen_max)])""",
    )
end

function (reducer::WsRvectorReducer)(operator::AbstractArray{T, 3}) where {T <: Number}
    size(operator, 3) == length(reducer.degens) ||
        throw(ArgumentError("Length of operator must match length of degens"))

    # Absorb R degeneracies into operator.
    operator_new = zeros(T, size(operator))
    for iR in axes(operator, 3)
        O = operator[:, :, iR]
        degen = reducer.degens[iR]
        operator_new[:, :, iR] = O ./ degen
    end
    return operator_new
end

function RvectorReducer(
        Rvectors::AbstractVector{<:AbstractVector{<:Integer}},
        Rdegens::AbstractVector{<:Integer},
    )
    return WsRvectorReducer(Rvectors, Rdegens)
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
    The number of tuples can be different for different `iR`.
    """
    mapping::Vector{Vector{Vec3{Int}}}

    """The combined R and T degeneracies of the original R vectors.
    `degens[:, :, iR0]` (`iR0 ∈ iR0s`) gives the degeneracy matrices of the original
    `iR0`-th R vector by combining the original `Rdegens` and the `Tdegens` of
    the T vectors.
    The `operator[:, :, iR0]` will be divided by `degens[:, :, iR0]` to absorb the degeneracies.
    """
    degens::Array{Int, 3}

    function MdrsRvectorReducer(
            Rvectors::AbstractVector, mapping::AbstractVector, degens::AbstractArray{<:Integer, 3}
        )
        length(Rvectors) == length(mapping) ||
            throw(ArgumentError("Length of Rvectors and mapping must match"))
        nR0 = size(degens, 3)
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
    return print(io, "MdrsRvectorReducer(n_Rvecs=", length(reducer.Rvectors), ")")
end

function Base.:(==)(a::MdrsRvectorReducer, b::MdrsRvectorReducer)
    return a.Rvectors == b.Rvectors && a.mapping == b.mapping && a.degens == b.degens
end

function Base.isequal(a::MdrsRvectorReducer, b::MdrsRvectorReducer)
    return (
        isequal(a.Rvectors, b.Rvectors) &&
            isequal(a.mapping, b.mapping) &&
            isequal(a.degens, b.degens)
    )
end

function Base.show(io::IO, ::MIME"text/plain", reducer::MdrsRvectorReducer)
    nR = length(reducer.Rvectors)
    nR0 = size(reducer.degens, 3)
    nw = (nR0 == 0) ? 0 : size(reducer.degens, 1)
    return print(
        io,
        """MdrsRvectorReducer(n_Rvecs = $(nR), n_Rvecs_original = $(nR0), n_wann = $(nw))""",
    )
end

function _validate_Tvectors_Tdegens(Tvectors::AbstractArray{<:Any, 3}, Tdegens::AbstractArray{<:Any, 3}, nR::Integer)
    nR == size(Tvectors, 3) == size(Tdegens, 3) ||
        throw(ArgumentError("number of R vectors mismatch between Tvectors and Tdegens"))
    nR > 0 || throw(ArgumentError("Length of Rvectors must be > 0"))
    # Number of Wannier functions
    nw = size(Tvectors, 1)
    nw == size(Tvectors, 2) == size(Tdegens, 1) == size(Tdegens, 2) ||
        throw(ArgumentError("n_wann mismatch between Tvectors and Tdegens"))
    for iR in 1:nR
        for Tvecs in Tvectors[:, :, iR]
            all(length.(Tvecs) .== 3) ||
                throw(ArgumentError("All T-vectors must be length 3"))
        end
    end
end

function MdrsRvectorReducer(
        Rvectors::AbstractVector{Vec3{IT}},
        Rdegens::AbstractVector{<:Integer},
        Tvectors::AbstractArray{<:AbstractVector{Vec3{IT}}, 3},
        Tdegens::AbstractArray{<:Integer, 3},
    ) where {IT <: Integer}
    # The length of original R-vectors
    nR0 = length(Rvectors)
    # The number of Wannier functions
    nw = size(Tvectors, 1)
    _validate_Rvectors_Rdegens(Rvectors, Rdegens)
    _validate_Tvectors_Tdegens(Tvectors, Tdegens, nR0)

    # Expanded R-vectors
    # 1. If empty, the new R+T will be added first.
    # Rvectors_new = Vector{Vec3{Int}}()
    # 2. Otherwise, we can first keep the same order as the orignal R-vectors in
    # the initial part, then the new R+T are added later.
    Rvectors_new = [Vec3{Int}(R) for R in Rvectors]

    # Mapping from new R vectors to old R-vector matrix entries.
    mapping = [Vec3{Int}[] for _ in eachindex(Rvectors_new)]
    # R and T degeneracies to absorb into operator.
    # Note the length is the same as the original R vectors, as the degeneracies
    # are applied to the original operator.
    degens = ones(Int, nw, nw, nR0)

    # generate expanded R vectors, which contains all the R+T
    for (iR0, R0) in enumerate(Rvectors)
        Tvecs = Tvectors[:, :, iR0]
        Nᴿ = Rdegens[iR0]
        # T degeneracy is a matrix
        Mᵀ = Tdegens[:, :, iR0]
        degens[:, :, iR0] .*= Nᴿ .* Mᵀ
        for n in axes(Mᵀ, 2)
            for m in axes(Mᵀ, 1)
                Nᵀ = Mᵀ[m, n]
                for iT in 1:Nᵀ
                    # The new R vector R = R0 + T
                    R = Vec3{Int}(R0 .+ Tvecs[m, n][iT])
                    iR = findfirst(==(R), Rvectors_new)
                    if isnothing(iR)
                        push!(Rvectors_new, R)
                        push!(mapping, [Vec3{Int}(iR0, m, n)])
                    else
                        push!(mapping[iR], Vec3{Int}(iR0, m, n))
                    end
                end
            end
        end
    end

    return MdrsRvectorReducer(Rvectors_new, mapping, degens)
end

function (reducer::MdrsRvectorReducer)(operator::AbstractArray{T, 3}) where {T <: Number}
    size(operator, 3) == size(reducer.degens, 3) ||
        throw(ArgumentError("Number of R-vectors mismatch between operator and degens"))

    # Simplified operator by absorbing R and T degeneracies
    # With length the same as the new R vectors
    nR = length(reducer.Rvectors)
    nw = size(operator, 1)
    operator_new = zeros(T, nw, nw, nR)
    for iR in 1:nR
        for (iR0, m, n) in reducer.mapping[iR]
            O = operator[:, :, iR0]
            degen = reducer.degens[:, :, iR0]
            operator_new[m, n, iR] += O[m, n] / degen[m, n]
        end
    end
    return operator_new
end

function RvectorReducer(
        Rvectors::AbstractVector{<:AbstractVector{<:Integer}},
        Rdegens::AbstractVector{<:Integer},
        Tvectors::AbstractArray{<:AbstractVector{<:AbstractVector{<:Integer}}, 3},
        Tdegens::AbstractArray{<:Integer, 3},
    )
    return MdrsRvectorReducer(Rvectors, Rdegens, Tvectors, Tdegens)
end
