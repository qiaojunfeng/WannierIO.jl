module WannierIOJLD2Ext

using WannierIO
using JLD2: JLD2

function WannierIO.write_w90_tb(
    filename::AbstractString,
    pack::WannierIO.OperatorPack,
    ::WannierIO.JLD2Format;
    atol::Real=0.0,
    index_type::Type{<:Integer}=Int32,
    value_type::Type{<:Number}=Float32,
    compress::Bool=true,
)
    opt = WannierIO.SparseOption(; atol, index_type, value_type)
    spack = WannierIO.sparsify(pack, opt)
    JLD2.jldsave(filename; pack=spack, compress)
    return nothing
end

function WannierIO.read_w90_tb(filename::AbstractString, ::WannierIO.JLD2Format)
    spack = JLD2.load(filename, "pack")
    return WannierIO.densify(spack)
end

end
