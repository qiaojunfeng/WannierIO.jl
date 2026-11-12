export CscPack, SparseOperatorPack

using SparseArrays: SparseMatrixCSC, sparse
using OrderedCollections: OrderedDict

"""
Pack a series of CSC matrices into a compact representation.

$(TYPEDEF)

This packs a dense `N × N × …` array, whose leading two axes are the matrix
axes, into one struct with concatenated nonzero values and indices. The axes
after the two matrix axes (components, R-vectors) are flattened in column-major
order into `L` matrices; each `N × N` matrix is stored in compressed sparse
column (CSC) format, and the `L` matrices are concatenated along the nonzero
dimension.

For example, `N = 3` and `L = 2` (`dims = [3, 3, 2]`):

```text
       l = 1             l = 2
    [ a  .  . ]       [ .  e  . ]
    [ .  .  b ]       [ .  .  . ]
    [ c  .  d ]       [ .  .  f ]

             |<---- l = 1 -->|<l = 2>|
  index      1   2   3   4   5   6   7
  nzval  = [ a   c   b   d   e   f ]
  rowval = [ 1   3   2   3   1   3 ]
  nzptr  = [ 1               5       7 ]

            l=1 l=2
  colptr = [ 1   1 ]    column j of matrix l is the local range
           [ 3   1 ]    colptr[j, l] : colptr[j+1, l] - 1, i.e. the
           [ 3   2 ]    indices nzptr[l] - 1 .+ (that range) of
           [ 5   3 ]    rowval and nzval
```

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

    """
    Dimensions of the dense array, `[N, N, trailing...]`, in column-major
    (Julia) order, e.g. `[n_wann, n_wann, n_Rvecs]` or
    `[n_wann, n_wann, 3, n_Rvecs]`. The product of the trailing dimensions
    is `L`.
    """
    dims::Vector{Int}

    function CscPack(
            nzptr::AbstractVector,
            colptr::AbstractMatrix,
            rowval::AbstractVector,
            nzval::AbstractVector,
            dims::AbstractVector{<:Integer},
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
        N = size(colptr, 1) - 1
        length(dims) >= 3 || error("invalid CSC pack: dims must have at least 3 entries")
        dims[1] == dims[2] == N || error("invalid CSC pack: dims $(dims) do not match N = $N")
        prod(dims[3:end]) == L ||
            error("invalid CSC pack: trailing dims $(dims[3:end]) do not multiply to L = $L")
        Ti = eltype(nzptr)
        Tv = eltype(nzval)
        return new{Tv, Ti}(
            nzptr, Matrix{Ti}(colptr), Vector{Ti}(rowval), nzval, Vector{Int}(dims)
        )
    end
end

"""
    CscPack{Tv, Ti}(mats, dims = [N, N, length(mats)])

Pack the sparse `N × N` matrices `mats`, the column-major flattening of a dense
array of dimensions `dims`.
"""
function CscPack{Tv, Ti}(
        mats::AbstractVector{<:SparseMatrixCSC},
        dims::AbstractVector{<:Integer} = _stacked_dims(mats),
    ) where {Tv <: Number, Ti <: Integer}
    L = length(mats)
    N = dims[1]
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
    return CscPack(nzptr, colptr, rowval, nzval, dims)
end

_stacked_dims(mats) = (N = isempty(mats) ? 0 : size(mats[1], 1); [N, N, length(mats)])

"""
    CscPack(op; atol = 0)
    CscPack{Tv}(op; atol = 0)
    CscPack{Tv, Ti}(op; atol = 0)

Pack the dense array `op`, `N × N × trailing...`, keeping the entries with
`abs(real(v)) > atol || abs(imag(v)) > atol`. In a kept entry, a real or
imaginary part of magnitude at most `atol` is set to zero. The default
`atol = 0` drops exactly the zeros, so the pack is lossless at value type
`Tv = eltype(op)`; the index type `Ti` defaults to `Int32`.

A real `Tv` for a complex `op` stores the real parts, and throws an
`InexactError` if an imaginary part larger than `atol` remains.

See `Array(pack)` for the inverse.
"""
CscPack(op::AbstractArray{T}; atol::Real = 0) where {T <: Number} = CscPack{T, Int32}(op; atol)

CscPack{Tv}(op::AbstractArray{<:Number}; atol::Real = 0) where {Tv <: Number} =
    CscPack{Tv, Int32}(op; atol)

function CscPack{Tv, Ti}(
        op::AbstractArray{<:Number}; atol::Real = 0
    ) where {Tv <: Number, Ti <: Integer}
    atol >= 0 || throw(ArgumentError("atol must be nonnegative, got $atol"))
    isempty(op) && throw(ArgumentError("operator must not be empty"))
    N = size(op, 1)
    mats = map(eachslice(_matrix_slices(op); dims = 3)) do O
        row = Ti[]
        col = Ti[]
        val = Tv[]
        @inbounds for n in 1:N, m in 1:N
            v = O[m, n]
            keep_re = abs(real(v)) > atol
            keep_im = abs(imag(v)) > atol
            (keep_re || keep_im) || continue
            vr = keep_re ? real(v) : zero(real(v))
            vi = keep_im ? imag(v) : zero(imag(v))
            push!(row, Ti(m))
            push!(col, Ti(n))
            push!(val, Tv(v isa Real ? vr : complex(vr, vi)))
        end
        return sparse(row, col, val, N, N)
    end
    return CscPack{Tv, Ti}(mats, collect(size(op)))
end

"""
    Array(pack::CscPack)
    Array{T}(pack::CscPack)

The dense array of `pack`, at its stored value type or converted to `T`.
"""
Base.Array(pack::CscPack{Tv}) where {Tv} = Array{Tv}(pack)

function Base.Array{T}(pack::CscPack) where {T <: Number}
    arr = zeros(T, size(pack))
    slices = _matrix_slices(arr)
    for l in 1:length(pack)
        slices[:, :, l] = pack[:, :, l]
    end
    return arr
end

# View an `N × N × trailing...` array as `N × N × L` with the trailing axes
# flattened in column-major order.
function _matrix_slices(arr::AbstractArray)
    ndims(arr) >= 3 || error("operator must have at least 3 dimensions, got $(ndims(arr))")
    size(arr, 1) == size(arr, 2) ||
        error("operator matrices must be square, got $(size(arr, 1)) × $(size(arr, 2))")
    return reshape(arr, size(arr, 1), size(arr, 2), :)
end

"Number of packed matrices `L`."
Base.length(p::CscPack) = size(p.colptr, 2)

"Size of the dense array, `(N, N, trailing...)`."
Base.size(p::CscPack) = Tuple(p.dims)

function Base.size(p::CscPack, d)
    return size(p)[d]
end

Base.ndims(p::CscPack) = length(p.dims)

"""
Return the `l`-th packed matrix as a `SparseMatrixCSC` with `p[:, :, l]`;
`l` is the column-major linear index over the trailing dimensions.
"""
function Base.getindex(
        p::CscPack{Tv, Ti}, ::Colon, ::Colon, l::Integer
    ) where {Tv <: Number, Ti <: Integer}
    1 <= l <= length(p) || throw(BoundsError(p, l))
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
          size: $(join(size(csc), "×"))
          total_nonzeros: $(nnz_total)
          avg_nonzeros_per_matrix: $(nnz_avg)
          value_type: $(eltype(csc.nzval))
          index_type: $(eltype(csc.nzptr))
        )""",
    )
end

# Whether every imaginary part of `arr` is at most `atol` in magnitude, so that
# `arr` can be stored as real.
_is_almost_real(arr::AbstractArray{<:Complex}, atol::Real) = all(v -> abs(imag(v)) <= atol, arr)
_is_almost_real(::AbstractArray{<:Real}, ::Real) = true

"""
Series of operators with sparse CSC storage.

# Fields
$(FIELDS)
"""
struct SparseOperatorPack{Tv <: Real, Ti <: Integer} <: AbstractOperatorPack
    "Short description of the operator set"
    header::String

    """Lattice vectors as columns of a 3×3 matrix, in Å. Always Float64,
    whatever the precision `Tv` of the operators."""
    lattice::Mat3{Float64}

    "R-vectors as columns of a `3 × n_Rvecs` matrix, integers."
    Rvectors::Matrix{Ti}

    "Mapping operator names to sparse packs."
    operators::OrderedDict{String, Union{CscPack{Tv, Ti}, CscPack{Complex{Tv}, Ti}}}

    function SparseOperatorPack(
            header::AbstractString,
            lattice::AbstractMatrix{<:Real},
            Rvectors::AbstractMatrix{Ti},
            operators::AbstractDict{<:AbstractString, <:CscPack},
        ) where {Ti <: Integer}
        Tv = _operator_real_type(operators)
        Rvs = Vec3{Ti}.(eachcol(Rvectors))
        _validate_operator_pack(Tv, Rvs, operators)

        Tops = Union{CscPack{Tv, Ti}, CscPack{Complex{Tv}, Ti}}
        ops = OrderedDict{String, Tops}()
        for (name, op) in pairs(operators)
            ops[String(name)] = op
        end

        return new{Tv, Ti}(String(header), Mat3{Float64}(lattice), Matrix{Ti}(Rvectors), ops)
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
    SparseOperatorPack(pack::OperatorPack; atol = 0)
    SparseOperatorPack{Tv}(pack::OperatorPack; atol = 0)
    SparseOperatorPack{Tv, Ti}(pack::OperatorPack; atol = 0)

Pack the operators of `pack` as [`CscPack`](@ref)s, for storage with
[`write_operator`](@ref).

- `Tv` is the real type of the stored values, by default that of `pack`, so
  that the default is lossless. `Tv = Float32` halves the values and rounds
  each element to about 7 significant digits.
- `Ti` is the integer type of the indices and R-vectors, by default `Int32`.
  It also counts the nonzeros of each operator, so a type smaller than `Int32`
  overflows for all but the smallest models.
- `atol` drops the entries whose real and imaginary parts are both at most
  `atol` in magnitude, and zeros either part when it alone is that small; see
  [`CscPack`](@ref). It is absolute, in the units of each operator (eV for a
  Hamiltonian, Å for a Berry connection). The default `0` drops only the exact
  zeros.

An operator whose imaginary parts are all at most `atol` is stored as real,
which at `atol = 0` loses nothing.

For a 70-Wannier-function Cu₂P₂Zr model with spin-orbit coupling (923
R-vectors, Hamiltonian and Berry connection), the HDF5 file is

| `Tv`      | `atol` | size      | max change in the bands |
|:----------|:-------|:----------|:------------------------|
| `Float64` | 0      | 110 MB    | 0                       |
| `Float32` | 0      | 72 MB     | 0.8 µeV                 |
| `Float32` | `1e-6` | 45 MB     | 16 µeV                  |

so `SparseOperatorPack{Float32}(pack; atol = 1e-6)` is a compact choice when
µeV accuracy is enough. See [`OperatorPack(::SparseOperatorPack)`](@ref) for
the inverse.
"""
SparseOperatorPack(pack::OperatorPack{Tv}; atol::Real = 0) where {Tv} =
    SparseOperatorPack{Tv, Int32}(pack; atol)

SparseOperatorPack{Tv}(pack::OperatorPack; atol::Real = 0) where {Tv <: Real} =
    SparseOperatorPack{Tv, Int32}(pack; atol)

function SparseOperatorPack{Tv, Ti}(
        pack::OperatorPack; atol::Real = 0
    ) where {Tv <: Real, Ti <: Integer}
    isempty(pack.operators) && throw(ArgumentError("operators must not be empty"))
    n_Rvectors(pack) > 0 || throw(ArgumentError("operators contain no R-vectors"))

    ops = OrderedDict(
        name => CscPack{_is_almost_real(op, atol) ? Tv : Complex{Tv}, Ti}(op; atol)
            for (name, op) in pairs(pack.operators)
    )
    Rvectors = Matrix{Ti}(reduce(hcat, pack.Rvectors))
    return SparseOperatorPack(pack.header, pack.lattice, Rvectors, ops)
end

"""
    OperatorPack(pack::SparseOperatorPack)
    OperatorPack{Tv}(pack::SparseOperatorPack)

The dense operators of `pack`, each at the value type it was stored with, or
with real type `Tv`. An operator stored as real comes back real.
"""
OperatorPack(pack::SparseOperatorPack{Tv}) where {Tv} = OperatorPack{Tv}(pack)

function OperatorPack{Tv}(pack::SparseOperatorPack{<:Real, Ti}) where {Tv <: Real, Ti}
    ops = OrderedDict(
        name => Array{eltype(op) <: Real ? Tv : Complex{Tv}}(op)
            for (name, op) in pairs(pack.operators)
    )
    Rvectors = Vec3{Ti}.(eachcol(pack.Rvectors))
    return OperatorPack(pack.header, pack.lattice, Rvectors, ops)
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
