module WannierIOHDF5Ext

using WannierIO
using HDF5: HDF5
using OrderedCollections: OrderedDict

const _DEFAULT_HDF5_DEFLATE = 4
const _DEFAULT_HDF5_SHUFFLE = true

function WannierIO.write_w90_tb(
    filename::AbstractString,
    pack::WannierIO.OperatorPack,
    fmt::WannierIO.HDF5Format;
    deflate::Integer=_DEFAULT_HDF5_DEFLATE,
    shuffle::Bool=_DEFAULT_HDF5_SHUFFLE,
    kwargs...,
)
    opt = WannierIO.SparseOption(; kwargs...)
    spack = WannierIO.sparsify(pack, opt)
    write_w90_tb(filename, spack, fmt; deflate, shuffle)
    return nothing
end

function WannierIO.write_w90_tb(
    filename::AbstractString,
    pack::WannierIO.SparseOperatorPack,
    ::WannierIO.HDF5Format;
    deflate::Integer=_DEFAULT_HDF5_DEFLATE,
    shuffle::Bool=_DEFAULT_HDF5_SHUFFLE,
)
    HDF5.h5open(filename, "w") do fid
        _h5write(fid, "header", pack.header; deflate, shuffle)
        _h5write(fid, "n_wann", pack.n_wann; deflate, shuffle)
        _h5write(fid, "n_Rvecs", pack.n_Rvecs; deflate, shuffle)
        _h5write(fid, "lattice", Matrix(pack.lattice); deflate, shuffle)
        _h5write(fid, "Rvectors", pack.Rvectors; deflate, shuffle)
        # Keep an explicit operator list for deterministic order and schema stability,
        # since HDF5 groups do not guarantee order.
        _h5write(fid, "operator_names", collect(keys(pack.operators)); deflate, shuffle)
        g = HDF5.create_group(fid, "operators")
        for (name, csc) in pack.operators
            _h5write_cscpack(HDF5.create_group(g, name), csc; deflate, shuffle)
        end
    end
    return nothing
end

function WannierIO.read_w90_tb(filename::AbstractString, ::WannierIO.HDF5Format)
    spack = HDF5.h5open(filename, "r") do fid
        header = read(fid["header"])
        n_wann = read(fid["n_wann"])
        n_Rvecs = read(fid["n_Rvecs"])
        lattice = read(fid["lattice"])
        Rvectors = read(fid["Rvectors"])
        op_names = read(fid["operator_names"])
        g = fid["operators"]
        operators = OrderedDict(name => _h5read_cscpack(g[name]) for name in op_names)

        Tr = eltype(lattice)
        p = WannierIO.SparseOperatorPack(
            header, WannierIO.Mat3{Tr}(lattice), Rvectors, operators
        )
        # Quick check
        p.n_wann == n_wann || error("n_wann mismatch: $n_wann vs $(p.n_wann)")
        p.n_Rvecs == n_Rvecs || error("n_Rvecs mismatch: $n_Rvecs vs $(p.n_Rvecs)")
        return p
    end
    return WannierIO.densify(spack)
end

function _h5write(fid, name::AbstractString, value; deflate::Integer, shuffle::Bool)
    if !(value isa AbstractArray)
        fid[name] = value
        return nothing
    end
    if eltype(value) <: AbstractString
        fid[name] = value
        return nothing
    end
    chunks = WannierIO._chunk_shape(value)
    fid[name, chunk = chunks, deflate = Int(deflate), shuffle = shuffle] = value
    return nothing
end

function _h5write_cscpack(
    g, csc::WannierIO.CscPack{V,I}; deflate::Integer, shuffle::Bool
) where {V,I}
    _h5write(g, "nzptr", csc.nzptr; deflate, shuffle)
    _h5write(g, "colptr", csc.colptr; deflate, shuffle)
    _h5write(g, "rowval", csc.rowval; deflate, shuffle)
    if V <: Complex
        _h5write(g, "nzval_re", real.(csc.nzval); deflate, shuffle)
        _h5write(g, "nzval_im", imag.(csc.nzval); deflate, shuffle)
    else
        _h5write(g, "nzval_re", csc.nzval; deflate, shuffle)
    end
    return nothing
end

function _h5read_cscpack(g)
    nzptr = read(g["nzptr"])
    colptr = read(g["colptr"])
    rowval = read(g["rowval"])
    nzval_re = read(g["nzval_re"])
    nzval = if haskey(g, "nzval_im")
        complex.(nzval_re, read(g["nzval_im"]))
    else
        nzval_re
    end
    return WannierIO.CscPack(nzptr, colptr, rowval, nzval)
end

end
