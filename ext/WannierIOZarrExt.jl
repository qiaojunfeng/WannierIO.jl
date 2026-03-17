module WannierIOZarrExt

using Mmap
using WannierIO
using Zarr: Zarr
using OrderedCollections: OrderedDict

const _DEFAULT_ZARR_CLEVEL = 5
const _DEFAULT_ZARR_SHUFFLE = true

function WannierIO.write_w90_tb(
    filename::AbstractString,
    pack::WannierIO.OperatorPack,
    fmt::Union{WannierIO.ZarrFormat,WannierIO.ZarrZipFormat};
    clevel::Integer=_DEFAULT_ZARR_CLEVEL,
    shuffle::Bool=_DEFAULT_ZARR_SHUFFLE,
    kwargs...,
)
    opt = WannierIO.SparseOption(; kwargs...)
    spack = WannierIO.sparsify(pack, opt)
    write_w90_tb(filename, spack, fmt; clevel, shuffle)
    return nothing
end

function WannierIO.write_w90_tb(
    filename::AbstractString,
    pack::WannierIO.SparseOperatorPack,
    ::WannierIO.ZarrFormat;
    clevel::Integer=_DEFAULT_ZARR_CLEVEL,
    shuffle::Bool=_DEFAULT_ZARR_SHUFFLE,
)
    ispath(filename) && error("File $filename already exists.")

    # I need to set the attrs at the creation of the group, otherwise Zarr.jl
    # does not write them.
    attrs = _zarr_tb_attrs(pack)
    # The default is a DirectoryStore, which creates a directory on disk.
    g = Zarr.zgroup(filename; attrs)

    _zarr_write_tb(g, pack; clevel, shuffle)

    return nothing
end

function WannierIO.write_w90_tb(
    filename::AbstractString,
    pack::WannierIO.SparseOperatorPack,
    ::WannierIO.ZarrZipFormat;
    clevel::Integer=_DEFAULT_ZARR_CLEVEL,
    shuffle::Bool=_DEFAULT_ZARR_SHUFFLE,
)
    ispath(filename) && error("File $filename already exists.")

    # I need to set the attrs at the creation of the group, otherwise Zarr.jl
    # does not write them.
    attrs = _zarr_tb_attrs(pack)
    # We want to store a single zip file, so we first use in-memory
    # DictStore and then zip it up at the end.
    store = Zarr.DictStore()
    g = Zarr.zgroup(store; attrs)

    _zarr_write_tb(g, pack; clevel, shuffle)

    open(filename, "w") do io
        Zarr.writezip(io, store)
    end
    return nothing
end

function _zarr_tb_attrs(pack::WannierIO.SparseOperatorPack)
    return Dict(
        "header" => pack.header,
        "n_wann" => pack.n_wann,
        "n_Rvecs" => pack.n_Rvecs,
        # Zarr does not preserve key/group order, so we store the
        # operator list as an explicit attribute for deterministic reading.
        "operator_names" => collect(keys(pack.operators)),
    )
end

function _zarr_write_tb(
    g::Zarr.ZGroup,
    pack::WannierIO.SparseOperatorPack;
    clevel::Integer=_DEFAULT_ZARR_CLEVEL,
    shuffle::Bool=_DEFAULT_ZARR_SHUFFLE,
)
    compressor = Zarr.BloscCompressor(; cname="zstd", clevel, shuffle)

    _zwrite(g, "lattice", Matrix(pack.lattice); compressor)
    _zwrite(g, "Rvectors", pack.Rvectors; compressor)

    ops_group = Zarr.zgroup(g, "operators")
    for (name, csc) in pack.operators
        _zwrite_cscpack(Zarr.zgroup(ops_group, name), csc; compressor)
    end

    return nothing
end

function WannierIO.read_w90_tb(filename::AbstractString, ::WannierIO.ZarrFormat)
    g = Zarr.zopen(filename)
    return _zarr_read_tb(g)
end

function WannierIO.read_w90_tb(filename::AbstractString, ::WannierIO.ZarrZipFormat)
    g = Zarr.zopen(Zarr.ZipStore(Mmap.mmap(filename)))
    return _zarr_read_tb(g)
end

function _zarr_read_tb(g::Zarr.ZGroup)
    header = g.attrs["header"]
    n_wann = g.attrs["n_wann"]
    n_Rvecs = g.attrs["n_Rvecs"]
    op_names = g.attrs["operator_names"]
    lattice = Array(g["lattice"][:, :])
    Rvectors = Array(g["Rvectors"][:, :])
    ops_group = g["operators"]
    operators = OrderedDict(name => _zread_cscpack(ops_group[name]) for name in op_names)

    Tr = eltype(lattice)
    spack = WannierIO.SparseOperatorPack(
        header, WannierIO.Mat3{Tr}(lattice), Rvectors, operators
    )
    # Just one quick sanity for dimensions
    spack.n_wann == n_wann || error("n_wann does not match operator dimensions")
    spack.n_Rvecs == n_Rvecs || error("n_Rvecs does not match operator dimensions")

    return WannierIO.densify(spack)
end

function _zwrite(g::Zarr.ZGroup, name::AbstractString, data::AbstractArray; compressor)
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

function _zwrite_cscpack(
    g::Zarr.ZGroup, csc::WannierIO.CscPack{V,I}; compressor
) where {V,I}
    _zwrite(g, "nzptr", csc.nzptr; compressor)
    _zwrite(g, "colptr", csc.colptr; compressor)
    _zwrite(g, "rowval", csc.rowval; compressor)
    # Requires Zarr.jl >= 0.9.5 for complex type, see
    # https://github.com/JuliaIO/Zarr.jl/pull/181
    _zwrite(g, "nzval", csc.nzval; compressor)
    return nothing
end

function _zread_cscpack(g::Zarr.ZGroup)
    nzptr = Array(g["nzptr"][:])
    colptr = Array(g["colptr"][:, :])
    rowval = Array(g["rowval"][:])
    nzval = Array(g["nzval"][:])
    return WannierIO.CscPack(nzptr, colptr, rowval, nzval)
end

end
