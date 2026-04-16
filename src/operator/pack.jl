export OperatorPack, pack, read_operator, write_operator

using OrderedCollections: OrderedDict

abstract type AbstractOperatorPack end

function n_wannier(operators::AbstractDict{String, <:Any})
    if isempty(operators)
        n_wann = 0
    else
        first_op = first(values(operators))
        n_wann = size(first_op, 1)
    end
    return n_wann
end

function _normalize_operators(::Type{Tv}, operators::AbstractDict) where {Tv <: Real}
    out = OrderedDict{String, Union{Array{Tv, 3}, Array{Complex{Tv}, 3}}}()
    for (name, op) in pairs(operators)
        key = String(name)
        if eltype(op) <: Real
            out[key] = Array{Tv, 3}(op)
        else
            out[key] = Array{Complex{Tv}, 3}(op)
        end
    end
    return out
end

function _validate_operator_pack(
        Tv::Type{<:Real},
        Rvectors::AbstractVector{<:Vec3{<:Integer}},
        operators::AbstractDict;
        n_Rvecs::Integer = length(Rvectors),
        n_wann::Union{Integer, Nothing} = nothing,
    )
    length(Rvectors) == n_Rvecs || error("length(Rvectors) != n_Rvecs")
    isnothing(n_wann) && (n_wann = n_wannier(operators))
    for (name, op) in pairs(operators)
        size(op, 3) == n_Rvecs ||
            error("length of operator `$name` ($(size(op, 3))) != n_Rvecs ($n_Rvecs)")
        To = eltype(op)
        To in (Tv, Complex{Tv}) ||
            error("operator `$name` has invalid element type $To, expected $Tv or $(Complex{Tv})")
        size(op)[1:2] == (n_wann, n_wann) ||
            error(
                "operator `$name` has invalid matrix size $(size(op)[1:2]), expected ($(n_wann), $(n_wann))",
            )
    end
    return nothing
end

"""
Series of operators with dense matrix storage.

The operators should have been simplified by [`WsRvectorReducer`](@ref)
or [`MdrsRvectorReducer`](@ref), i.e., no `Rdegens` or `Tvectors`.

# Fields
$(FIELDS)
"""
struct OperatorPack{Tv <: Real, Ti <: Integer} <: AbstractOperatorPack
    "Short description of the operator set"
    header::String

    "Lattice vectors as columns of a 3×3 matrix, in Å."
    lattice::Mat3{Tv}

    "R-vectors as a vector of 3-component integer vectors."
    Rvectors::Vector{Vec3{Ti}}

    """Mapping operator names to vectors of dense matrices.
    The operators can be either real or complex.
    """
    operators::OrderedDict{String, Union{Array{Tv, 3}, Array{Complex{Tv}, 3}}}

    function OperatorPack(
            header::AbstractString,
            lattice::AbstractMatrix{Tv},
            Rvectors::AbstractVector{<:Vec3{Ti}},
            operators::AbstractDict,
        ) where {Tv <: Real, Ti <: Integer}
        header = String(header)
        lattice = Mat3{Tv}(lattice)
        _validate_operator_pack(Tv, Rvectors, operators)
        Rvectors = collect(Vec3{Ti}.(Rvectors))
        operators = _normalize_operators(Tv, operators)
        return new{Tv, Ti}(header, lattice, Rvectors, operators)
    end
end

n_Rvectors(op::OperatorPack) = length(op.Rvectors)
n_wannier(op::OperatorPack) = n_wannier(op.operators)

function Base.show(io::IO, op::OperatorPack)
    n_ops = length(op.operators)
    return print(io, "OperatorPack(n_Rvecs=$(n_Rvectors(op)), n_wann=$(n_wannier(op)), n_operators=$(n_ops))")
end

function Base.show(io::IO, ::MIME"text/plain", op::OperatorPack)
    n_ops = length(op.operators)
    op_names = collect(keys(op.operators))

    return print(
        io,
        """OperatorPack(
          header: $(op.header)
          n_Rvecs: $(n_Rvectors(op))
          n_wann: $(n_wannier(op))
          operators ($(n_ops)): $(join(op_names, ", "))
        )""",
    )
end

"""
Pack tight-binding operators into an [`OperatorPack`](@ref) struct.

e.g., from [`TbDat`](@ref) and [`WsvecDat`](@ref) structs to [`OperatorPack`](@ref).
"""
function pack end

"""
    read_operator(file, format)
    read_operator(file)

Read operators from a backend storage format (HDF5/JLD2/Zarr) and return an
[`OperatorPack`](@ref).

When called without a format argument, the format is inferred from the file
extension via [`detect_operator_format`](@ref).
"""
function read_operator end

"""
    write_operator(file, pack, format; kwargs...)
    write_operator(file, pack; kwargs...)

Write an operator pack to a backend storage format (HDF5/JLD2/Zarr).

When called without a format argument, the format is inferred from the file
extension via [`detect_operator_format`](@ref).
"""
function write_operator end

function read_operator(file::AbstractString, fmt::AbstractFileFormat)
    error(
        "Format `$(format_name(fmt))` requires loading the corresponding package. " *
            "See the WannierIO documentation for supported formats.",
    )
end

function read_operator(file::AbstractString)
    return read_operator(file, detect_operator_format(file))
end

function write_operator(
        file::AbstractString, pack::AbstractOperatorPack, fmt::AbstractFileFormat; kwargs...
    )
    error(
        "Format `$(format_name(fmt))` requires loading the corresponding package. " *
            "See the WannierIO documentation for supported formats.",
    )
end

function write_operator(file::AbstractString, pack::AbstractOperatorPack; kwargs...)
    return write_operator(file, pack, detect_operator_format(file); kwargs...)
end
