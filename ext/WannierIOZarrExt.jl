module WannierIOZarrExt

using WannierIO
using Zarr: Zarr
using OrderedCollections: OrderedDict

function WannierIO.write_w90_tb(
    path::AbstractString,
    pack::WannierIO.OperatorPack,
    ::WannierIO.ZarrFormat;
    atol::Real=0.0,
    index_type::Type{<:Integer}=Int32,
    value_type::Type{<:Number}=Float32,
    clevel::Integer=5,
)
    opt = WannierIO.SparseOption(; atol, index_type, value_type)
    spack = WannierIO.sparsify(pack, opt)
    compressor = Zarr.BloscCompressor(0, Int(clevel), "zstd", true)

    g = Zarr.zgroup(path)
    _zwrite(g, "header", [spack.header]; compressor)
    _zwrite(g, "n_wann", reshape([spack.n_wann], 1); compressor)
    _zwrite(g, "n_Rvecs", reshape([spack.n_Rvecs], 1); compressor)
    _zwrite(g, "lattice", Matrix(spack.lattice); compressor)
    _zwrite(g, "Rvectors", spack.Rvectors; compressor)
    _zwrite(g, "operator_names", collect(keys(spack.operators)); compressor)

    ops_group = Zarr.zgroup(g, "operators")
    for (name, csc) in spack.operators
        _zwrite_cscpack(Zarr.zgroup(ops_group, name), csc; compressor)
    end

    return nothing
end

function WannierIO.read_w90_tb(path::AbstractString, ::WannierIO.ZarrFormat)
    g = Zarr.zopen(path)

    header = String(only(g["header"][:]))
    n_wann = only(g["n_wann"][:])
    n_Rvecs = only(g["n_Rvecs"][:])
    lattice = Array(g["lattice"][:, :])
    Rvectors = Array(g["Rvectors"][:, :])
    op_names = String.(g["operator_names"][:])
    ops_group = g["operators"]
    operators = OrderedDict(name => _zread_cscpack(ops_group[name]) for name in op_names)

    Tr = eltype(lattice)
    spack = WannierIO.SparseOperatorPack(
        header, WannierIO.Mat3{Tr}(lattice), Rvectors, operators
    )
    return WannierIO.densify(spack)
end

function _zwrite(g, name::AbstractString, data::AbstractArray; compressor)
    chunks = WannierIO._chunk_shape(data)
    z = Zarr.zcreate(eltype(data), g, name, size(data)...; chunks, compressor)
    if ndims(data) == 1
        z[:] = data
    elseif ndims(data) == 2
        z[:, :] = data
    else
        error("unsupported ndims for zarr writer")
    end
    return nothing
end

function _zwrite_cscpack(g, csc::WannierIO.CscPack{V,I}; compressor) where {V,I}
    _zwrite(g, "nzptr", csc.nzptr; compressor)
    _zwrite(g, "colptr", csc.colptr; compressor)
    _zwrite(g, "rowval", csc.rowval; compressor)
    if V <: Complex
        _zwrite(g, "nzval_re", real.(csc.nzval); compressor)
        _zwrite(g, "nzval_im", imag.(csc.nzval); compressor)
    else
        _zwrite(g, "nzval_re", csc.nzval; compressor)
    end
    return nothing
end

function _zread_cscpack(g)
    nzptr = Array(g["nzptr"][:])
    colptr = Array(g["colptr"][:, :])
    rowval = Array(g["rowval"][:])
    nzval_re = Array(g["nzval_re"][:])
    nzval = if haskey(g, "nzval_im")
        complex.(nzval_re, Array(g["nzval_im"][:]))
    else
        nzval_re
    end
    return WannierIO.CscPack(nzptr, colptr, rowval, nzval)
end

end
