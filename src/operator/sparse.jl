export sparsify, densify

using SparseArrays: SparseMatrixCSC, sparse
using OrderedCollections: OrderedDict

"""
Pack a series of CSC matrices into a compact representation.

$(TYPEDEF)

This packs a `N × N × L` sparse array into one struct with concatenated
nonzero values and indices. Each `N × N` matrix is stored in compressed sparse
column (CSC) format, and the `L` matrices are concatenated along the nonzero dimension.

# Fields

$(FIELDS)
"""
struct CscPack{Tv <: Number, Ti <: Integer}
    """
    1-based offsets into concatenated nonzero buffers for each matrix.
    Length is `L + 1`, and `nzptr[l]:nzptr[l+1]-1` are the entries of
    the `l`-th matrix in `rowval` and `nzval`.
    """
    nzptr::Vector{Ti}

    """
    CSC column pointer table for each packed matrix.
    Size is `(N + 1, L)`, and `colptr[:, l]` is the column pointer of
    the `l`-th matrix.
    """
    colptr::Matrix{Ti}

    """
    Concatenated row indices of all nonzeros.
    Length is the total number of nonzeros across all matrices.
    """
    rowval::Vector{Ti}

    """
    Concatenated nonzero values of all matrices.
    Length is the total number of nonzeros across all matrices.
    """
    nzval::Vector{Tv}

    function CscPack(
            nzptr::AbstractVector,
            colptr::AbstractMatrix,
            rowval::AbstractVector,
            nzval::AbstractVector,
        )
        L = length(nzptr) - 1
        L >= 0 || error("invalid CSC pack: nzptr must have length >= 1")
        size(colptr, 2) == L || error("invalid CSC pack: colptr length mismatch")
        length(rowval) == length(nzval) || error("invalid CSC pack: rowval/nzval mismatch")

        nzptr[1] == 1 || error("invalid CSC pack: nzptr must start at 1")
        nzptr[end] == length(rowval) + 1 || error("invalid CSC pack: nzptr end mismatch")
        @inbounds for l in 1:L
            nzptr[l] <= nzptr[l + 1] ||
                error("invalid CSC pack: nzptr must be nondecreasing")
        end

        Ti = eltype(nzptr)
        Tv = eltype(nzval)
        return new{Tv, Ti}(nzptr, Matrix{Ti}(colptr), Vector{Ti}(rowval), nzval)
    end
end

function CscPack{Tv, Ti}(
        mats::AbstractVector{<:SparseMatrixCSC}
    ) where {Tv <: Number, Ti <: Integer}
    L = length(mats)
    N = isempty(mats) ? 0 : size(mats[1], 1)

    nzptr = Vector{Ti}(undef, L + 1)
    nzptr[1] = 1
    colptr = Matrix{Ti}(undef, N + 1, L)
    rowval = Vector{Ti}()
    nzval = Vector{Tv}()

    @inbounds for l in 1:L
        M = mats[l]
        size(M, 1) == N || error("inconsistent matrix shape")
        size(M, 2) == N || error("inconsistent matrix shape")

        colptr[:, l] = Ti.(M.colptr)
        append!(rowval, Ti.(M.rowval))
        append!(nzval, Tv.(M.nzval))
        nzptr[l + 1] = Ti(length(rowval) + 1)
    end

    return CscPack(nzptr, colptr, rowval, nzval)
end

function CscPack{Tv, Ti}(arr::AbstractArray{T, 3}) where {Tv <: Number, Ti <: Integer, T <: Number}
    mats = [sparse(arr[:, :, i]) for i in axes(arr, 3)]
    return CscPack{Tv, Ti}(mats)
end

"Number of packed matrices."
Base.length(p::CscPack) = size(p.colptr, 2)

"Size `N × N × L`."
function Base.size(p::CscPack)
    N = size(p.colptr, 1) - 1
    L = size(p.colptr, 2)
    return (N, N, L)
end

function Base.size(p::CscPack, d)
    size(p)[d]
end

"Return the `l`-th matrix as a `SparseMatrixCSC` with `p[:, :, l]`."
function Base.getindex(
        p::CscPack{Tv, Ti}, ::Colon, ::Colon, l::Integer
    ) where {Tv <: Number, Ti <: Integer}
    1 <= l <= size(p, 3) || throw(BoundsError(p, l))
    r = p.nzptr[l]:(p.nzptr[l + 1] - 1)
    N = size(p, 1)
    return SparseMatrixCSC{Tv, Ti}(N, N, p.colptr[:, l], p.rowval[r], p.nzval[r])
end

Base.eltype(::Type{CscPack{Tv, Ti}}) where {Tv <: Number, Ti <: Integer} = Tv
Base.isempty(p::CscPack) = prod(size(p)) == 0

function Base.show(io::IO, csc::CscPack)
    L = length(csc)
    N = size(csc, 1)
    nnz_total = length(csc.nzval)
    return print(io, "CscPack(L=$(L), N=$(N), nnz=$(nnz_total))")
end

function Base.show(io::IO, ::MIME"text/plain", csc::CscPack)
    L = length(csc)
    N = size(csc, 1)
    nnz_total = length(csc.nzval)
    nnz_avg = L == 0 ? 0 : div(nnz_total, L)

    return print(
        io,
        """CscPack(
          size: $(N)×$(N)×$(L)
          total_nonzeros: $(nnz_total)
          avg_nonzeros_per_matrix: $(nnz_avg)
          value_type: $(eltype(csc.nzval))
          index_type: $(eltype(csc.nzptr))
        )""",
    )
end

"""
Control sparsification behavior.

# Fields
$(FIELDS)
"""
Base.@kwdef struct SparseOption
    """Absolute tolerance for treating values as zero. Must be nonnegative.
    Elements with `abs(real(v)) < atol` and `abs(imag(v)) < atol` are dropped.
    """
    atol::Float64 = 1.0e-6

    "Integer type for storing row and column indices in sparse matrices."
    index_type::Type{<:Integer} = Int32

    """Numeric type for storing nonzero values in sparse matrices.
    Can be real or complex. If real, will attempt to convert complex values
    to real if all imaginary parts are smaller than `atol`."""
    value_type::Type{<:Number} = Float32

    function SparseOption(atol::Real, index_type, value_type)
        atol >= 0 || error("atol must be nonnegative")
        return new(Float64(atol), index_type, value_type)
    end
end

"""
Check if the imaginary part is negligible compared to `atol`.
"""
function is_almost_real(arr::AbstractArray{T, 3}, atol::Real) where {T <: Complex}
    @inbounds for M in arr
        for v in M
            abs(imag(v)) < atol || return false
        end
    end
    return true
end

is_almost_real(::AbstractArray{<:Real, 3}, ::Real) = true

"""
    _force_value_type(operator_eltype, value_type)

Return the sparse storage element type for an operator with element type
`operator_eltype` and requested precision `value_type`.
"""
_force_value_type(::Type{T}, ::Type{S}) where {T <: Real, S <: Number} = S
function _force_value_type(::Type{Complex{T}}, ::Type{S}) where {T <: Real, S <: Number}
    return Complex{real(S)}
end

"""
    $(SIGNATURES)

Convert an operator to sparse CSC representation.

Only entries with `abs(real(v)) >= atol || abs(imag(v)) >= atol` are stored.

See [`densify`](@ref) for the inverse.

# Arguments
- `op`: a `n_wann × n_wann × n_Rvecs` array representing the operator at different R-vectors.

# Keyword Arguments
"""
function sparsify(op::AbstractArray{T, 3}, opt::SparseOption = SparseOption()) where {T <: Number}
    isempty(op) && error("operator list must not be empty")

    Ti = opt.index_type
    N = size(op, 1)
    to_real = (T <: Complex) && (opt.value_type <: Real) && is_almost_real(op, opt.atol)
    Tvout = to_real ? opt.value_type : _force_value_type(T, opt.value_type)
    Trout = real(Tvout)

    mats = map(1:size(op, 3)) do i
        O = op[:, :, i]
        size(O) == (N, N) || error("operator matrices must have size ($N, $N)")
        row = Ti[]
        col = Ti[]
        val = Tvout[]
        @inbounds for n in 1:N
            for m in 1:N
                v = O[m, n]
                rv = real(v)
                iv = imag(v)
                if abs(rv) >= opt.atol || abs(iv) >= opt.atol
                    push!(row, Ti(m))
                    push!(col, Ti(n))
                    ro = abs(rv) >= opt.atol ? Trout(rv) : Trout(0)
                    io = abs(iv) >= opt.atol ? Trout(iv) : Trout(0)
                    push!(val, to_real ? ro : Complex(ro, io))
                end
            end
        end
        return sparse(row, col, val, N, N)
    end

    return CscPack{Tvout, Ti}(mats)
end

"""
    $(SIGNATURES)

Reconstruct dense matrices from sparse CSC matrices produced by [`sparsify`](@ref).

Set `value_type` to force conversion. Use `nothing` to preserve stored type.
"""
function densify(pack::CscPack{Tv, Ti}; value_type::Type{<:Number} = Tv) where {Tv, Ti}
    L = size(pack, 3)
    T = _force_value_type(Tv, value_type)
    arr = Array{T, 3}(undef, size(pack))
    for l in 1:L
        arr[:, :, l] = pack[:, :, l]
    end
    return arr
end

"""
Series of operators with sparse CSC storage.

# Fields
$(FIELDS)
"""
struct SparseOperatorPack{Tv <: Real, Ti <: Integer} <: AbstractOperatorPack
    "Short description of the operator set"
    header::String

    "Lattice vectors as columns of a 3×3 matrix, in Å."
    lattice::Mat3{Tv}

    "R-vectors as columns of a `3 × n_Rvecs` matrix, integers."
    Rvectors::Matrix{Ti}

    "Mapping operator names to sparse packs."
    operators::OrderedDict{String, Union{CscPack{Tv, Ti}, CscPack{Complex{Tv}, Ti}}}

    function SparseOperatorPack(
            header::AbstractString,
            lattice::AbstractMatrix{Tv},
            Rvectors::AbstractMatrix{Ti},
            operators::AbstractDict{<:AbstractString, <:CscPack},
        ) where {Tv <: Real, Ti <: Integer}
        Rvs = Vec3{Ti}.(eachcol(Rvectors))
        _validate_operator_pack(Tv, Rvs, operators)

        Tops = Union{CscPack{Tv, Ti}, CscPack{Complex{Tv}, Ti}}
        ops = OrderedDict{String, Tops}()
        for (name, op) in pairs(operators)
            ops[String(name)] = op
        end

        return new{Tv, Ti}(String(header), Mat3{Tv}(lattice), Matrix{Ti}(Rvectors), ops)
    end
end

n_Rvectors(pack::SparseOperatorPack) = size(pack.Rvectors, 2)
n_wannier(pack::SparseOperatorPack) = n_wannier(pack.operators)

function Base.show(io::IO, sop::SparseOperatorPack)
    n_ops = length(sop.operators)
    return print(
        io,
        "SparseOperatorPack(n_Rvecs=$(n_Rvectors(sop)), n_wann=$(n_wannier(sop)), n_operators=$(n_ops))",
    )
end

function Base.show(io::IO, ::MIME"text/plain", sop::SparseOperatorPack)
    n_ops = length(sop.operators)
    op_names = collect(keys(sop.operators))

    return print(
        io,
        """SparseOperatorPack(
          header: $(sop.header)
          n_Rvecs: $(n_Rvectors(sop))
          n_wann: $(n_wannier(sop))
          operators ($(n_ops)): $(join(op_names, ", "))
        )""",
    )
end

"""
    $(SIGNATURES)

Pack a set of Wannier90 operators into sparse CSC matrices suitable for
storage in HDF5, JLD2, or Zarr.

See [`densify`](@ref) for the inverse.
"""
function sparsify(pack::OperatorPack, opt::SparseOption = SparseOption())
    isempty(pack.operators) && error("operators must not be empty")
    n_Rvectors(pack) > 0 || error("operators contain no R-vectors")

    ops = OrderedDict(name => sparsify(op, opt) for (name, op) in pairs(pack.operators))

    Tr = real(opt.value_type)
    Ti = opt.index_type
    Rvectors = if isempty(pack.Rvectors)
        Matrix{Ti}(undef, 3, 0)
    else
        Matrix{Ti}(reduce(hcat, Vec3{Ti}.(pack.Rvectors)))
    end
    return SparseOperatorPack(pack.header, Mat3{Tr}(pack.lattice), Rvectors, ops)
end

"""
    $(SIGNATURES)

Reconstruct dense operators from a sparse CSC pack produced by
[`sparsify`](@ref).
"""
function densify(
        pack::SparseOperatorPack{Tr, I}; value_type::Type{<:Number} = Complex{Tr}
    ) where {Tr, I}
    ops = OrderedDict(
        name => densify(op; value_type) for (name, op) in pairs(pack.operators)
    )

    Trout = real(value_type)
    Rvectors = Vec3{I}.(eachcol(pack.Rvectors))
    return OperatorPack(pack.header, Mat3{Trout}(pack.lattice), Rvectors, ops)
end

# 2^18 = 262_144 elements max per 1D chunk.
const _MAX_CHUNK_ELEMS_1D = 1 << 18
# 2^9 = 512 per edge; 512^2 matches the 1D element budget.
const _MAX_CHUNK_EDGE_2D = 1 << 9

"""
    _chunk_shape(data)

Return a conservative HDF5/Zarr chunk shape.

For 1D arrays, cap the chunk at `_MAX_CHUNK_ELEMS_1D` elements.
For 2D arrays, cap each edge at `_MAX_CHUNK_EDGE_2D`, i.e. at most
`_MAX_CHUNK_EDGE_2D^2 == _MAX_CHUNK_ELEMS_1D` elements per chunk.
"""
function _chunk_shape(data::AbstractArray)
    if ndims(data) == 1
        return (max(1, min(length(data), _MAX_CHUNK_ELEMS_1D)),)
    elseif ndims(data) == 2
        m, n = size(data)
        return (max(1, min(m, _MAX_CHUNK_EDGE_2D)), max(1, min(n, _MAX_CHUNK_EDGE_2D)))
    else
        return size(data)
    end
end
