export OperatorPack, pack

using OrderedCollections: OrderedDict

function _infer_n_wann(operators::AbstractDict)
    if isempty(operators)
        n_wann = 0
    else
        first_op = first(values(operators))
        if isempty(first_op)
            n_wann = 0
        else
            O = first(first_op)
            size(O, 1) == size(O, 2) || error("operator matrices must be square")
            n_wann = size(O, 1)
        end
    end
    return n_wann
end

function _normalize_operators(::Type{Tv}, operators::AbstractDict) where {Tv<:Real}
    Tval = Union{Vector{Matrix{Tv}},Vector{Matrix{Complex{Tv}}}}
    out = OrderedDict{String,Tval}()
    for (name, op) in pairs(operators)
        key = String(name)
        if eltype(eltype(op)) <: Real
            out[key] = [Matrix{Tv}(O) for O in op]
        else
            out[key] = [Matrix{Complex{Tv}}(O) for O in op]
        end
    end
    return out
end

function _validate_operator_pack(
    Tv::Type{<:Real},
    n_wann::Integer,
    n_Rvecs::Integer,
    Rvectors::AbstractVector{<:Vec3{<:Integer}},
    operators::AbstractDict,
)
    length(Rvectors) == n_Rvecs || error("length(Rvectors) != n_Rvecs")
    for (name, op) in pairs(operators)
        length(op) == n_Rvecs ||
            error("length of operator `$name` ($(length(op))) != n_Rvecs ($n_Rvecs)")
        To = eltype(eltype(op))
        To in (Tv, Complex{Tv}) ||
            error("operator `$name` has invalid element type $To, expected $Tv or $Tc")
        for O in op
            size(O) == (n_wann, n_wann) || error(
                "operator `$name` has invalid matrix size $(size(O)), expected ($(n_wann), $(n_wann))",
            )
        end
    end
end

"""
Series of operators with dense matrix storage.

The operators should have been simplified by [`WsRvectorReducer`](@ref)
or [`MdrsRvectorReducer`](@ref), i.e., no `Rdegens` or `Tvectors`.

# Fields
$(FIELDS)
"""
struct OperatorPack{Tv<:Real,Ti<:Integer}
    "Short description of the operator set"
    header::String

    "Lattice vectors as columns of a 3×3 matrix, in Å."
    lattice::Mat3{Tv}

    "R-vectors as a vector of 3-component integer vectors."
    Rvectors::Vector{Vec3{Ti}}

    """Mapping operator names to vectors of dense matrices.
    The operators can be either real or complex.
    """
    operators::OrderedDict{String,Union{Vector{Matrix{Tv}},Vector{Matrix{Complex{Tv}}}}}

    "Number of R-vectors."
    n_Rvecs::Ti

    "Number of Wannier functions."
    n_wann::Ti

    function OperatorPack(
        header::AbstractString,
        lattice::AbstractMatrix{Tv},
        Rvectors::AbstractVector{<:Vec3{Ti}},
        operators::AbstractDict,
    ) where {Tv<:Real,Ti<:Integer}
        header = String(header)
        lattice = Mat3{Tv}(lattice)
        n_Rvecs = Ti(length(Rvectors))
        n_wann = Ti(_infer_n_wann(operators))
        _validate_operator_pack(Tv, n_wann, n_Rvecs, Rvectors, operators)
        Rvectors = collect(Vec3{Ti}.(Rvectors))
        operators = _normalize_operators(Tv, operators)
        return new{Tv,Ti}(header, lattice, Rvectors, operators, n_Rvecs, n_wann)
    end
end

function Base.show(io::IO, op::OperatorPack)
    n_ops = length(op.operators)
    print(
        io, "OperatorPack(n_Rvecs=$(op.n_Rvecs), n_wann=$(op.n_wann), n_operators=$(n_ops))"
    )
end

function Base.show(io::IO, ::MIME"text/plain", op::OperatorPack)
    n_ops = length(op.operators)
    op_names = collect(keys(op.operators))

    print(
        io,
        """OperatorPack(
          header: $(op.header)
          n_Rvecs: $(op.n_Rvecs)
          n_wann: $(op.n_wann)
          operators ($(n_ops)): $(join(op_names, ", "))
        )""",
    )
end

"""
Pack tight-binding operators into an [`OperatorPack`](@ref) struct.

e.g., from [`TbDat`](@ref) and [`WsvecData`](@ref) structs to [`OperatorPack`](@ref).
"""
function pack end
