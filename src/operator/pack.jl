export OperatorPack, read_operator, write_operator

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
    out = OrderedDict{String, Union{Array{Tv}, Array{Complex{Tv}}}}()
    for (name, op) in pairs(operators)
        key = String(name)
        if eltype(op) <: Real
            out[key] = Array{Tv}(op)
        else
            out[key] = Array{Complex{Tv}}(op)
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
        dims = size(op)
        length(dims) >= 3 || error(
            "operator `$name` must be n_wann × n_wann × … × n_Rvecs, got size $dims",
        )
        dims[end] == n_Rvecs ||
            error("length of operator `$name` ($(dims[end])) != n_Rvecs ($n_Rvecs)")
        To = eltype(op)
        To in (Tv, Complex{Tv}) ||
            error("operator `$name` has invalid element type $To, expected $Tv or $(Complex{Tv})")
        dims[1:2] == (n_wann, n_wann) ||
            error(
            "operator `$name` has invalid matrix size $(dims[1:2]), expected ($(n_wann), $(n_wann))",
        )
    end
    return nothing
end

"""
Named Wannier operators on a common set of R-vectors, with dense storage.

An `OperatorPack` is a storage container: it checks shapes, but holds no
physics. Each operator is an array `O[m, n, components..., iR]` of
`n_wann × n_wann × component dimensions × n_Rvecs`, e.g. a Hamiltonian
`n_wann × n_wann × n_Rvecs` or a vector operator
`n_wann × n_wann × 3 × n_Rvecs`. The R-vectors are unique and carry no
degeneracy weights, so an operator at a kpoint is the plain sum
``O(\\mathbf{k}) = \\sum_{\\mathbf{R}} e^{i 2\\pi \\mathbf{k} \\cdot \\mathbf{R}} O(\\mathbf{R})``
(`k` and `R` fractional). Wigner--Seitz degeneracies and minimum-distance
replicas of Wannier90 `tb.dat`/`wsvec.dat` files must have been absorbed
before packing, e.g. by Wannier.jl's `OperatorPack(::TightBindingModel)`.

The operator names are free-form strings; the dense layout is preserved by
[`sparsify`](@ref)/[`densify`](@ref) and by [`write_operator`](@ref)/
[`read_operator`](@ref).

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

    """Mapping operator names to dense `n_wann × n_wann × … × n_Rvecs` arrays.
    The operators can be either real or complex.
    """
    operators::OrderedDict{String, Union{Array{Tv}, Array{Complex{Tv}}}}

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
