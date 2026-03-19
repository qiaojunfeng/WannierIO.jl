module WannierIOJLD2Ext

using WannierIO
using JLD2: JLD2

const _DEFAULT_JLD2_COMPRESS = true

function WannierIO.write_operators(
    filename::AbstractString,
    pack::WannierIO.OperatorPack,
    fmt::WannierIO.JLD2Format;
    compress::Bool=_DEFAULT_JLD2_COMPRESS,
    kwargs...,
)
    opt = WannierIO.SparseOption(; kwargs...)
    spack = WannierIO.sparsify(pack, opt)
    write_operators(filename, spack, fmt; compress)
    return nothing
end

function WannierIO.write_operators(
    filename::AbstractString,
    pack::WannierIO.SparseOperatorPack,
    ::WannierIO.JLD2Format;
    compress::Bool=_DEFAULT_JLD2_COMPRESS,
)
    JLD2.jldsave(filename, compress; pack)
end

function WannierIO.read_operators(filename::AbstractString, ::WannierIO.JLD2Format)
    spack = JLD2.load(filename, "pack")
    return WannierIO.densify(spack)
end

end
