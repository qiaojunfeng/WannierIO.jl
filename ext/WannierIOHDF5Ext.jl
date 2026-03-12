module WannierIOHDF5Ext

using WannierIO
using HDF5: HDF5
using OrderedCollections: OrderedDict

function WannierIO.write_w90_tb(
    filename::AbstractString,
    pack::WannierIO.OperatorPack,
    ::WannierIO.HDF5Format;
    atol::Real=0.0,
    index_type::Type{<:Integer}=Int32,
    value_type::Type{<:Number}=Float32,
    deflate::Integer=4,
)
    opt = WannierIO.SparseOption(; atol, index_type, value_type)
    spack = WannierIO.sparsify(pack, opt)
    HDF5.h5open(filename, "w") do fid
        _h5write(fid, "header", spack.header; deflate)
        _h5write(fid, "n_wann", spack.n_wann; deflate)
        _h5write(fid, "n_Rvecs", spack.n_Rvecs; deflate)
        _h5write(fid, "lattice", Matrix(spack.lattice); deflate)
        _h5write(fid, "Rvectors", spack.Rvectors; deflate)
        _h5write(fid, "operator_names", collect(keys(spack.operators)); deflate)
        g = HDF5.create_group(fid, "operators")
        for (name, csc) in spack.operators
            _h5write_cscpack(HDF5.create_group(g, name), csc; deflate)
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
        WannierIO.SparseOperatorPack(
            header, WannierIO.Mat3{Tr}(lattice), Rvectors, operators
        )
    end
    return WannierIO.densify(spack)
end

function _h5write(fid, name::AbstractString, value; deflate::Integer)
    if !(value isa AbstractArray)
        fid[name] = value
        return nothing
    end
    if eltype(value) <: AbstractString
        fid[name] = value
        return nothing
    end
    chunks = WannierIO._chunk_shape(value)
    fid[name, chunk = chunks, deflate = Int(deflate)] = value
    return nothing
end

function _h5write_cscpack(g, csc::WannierIO.CscPack{V,I}; deflate::Integer) where {V,I}
    _h5write(g, "nzptr", csc.nzptr; deflate)
    _h5write(g, "colptr", csc.colptr; deflate)
    _h5write(g, "rowval", csc.rowval; deflate)
    if V <: Complex
        _h5write(g, "nzval_re", real.(csc.nzval); deflate)
        _h5write(g, "nzval_im", imag.(csc.nzval); deflate)
    else
        _h5write(g, "nzval_re", csc.nzval; deflate)
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
