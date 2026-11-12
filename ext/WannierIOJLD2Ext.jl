module WannierIOJLD2Ext

using WannierIO
using JLD2: JLD2

const _DEFAULT_JLD2_COMPRESS = true

function WannierIO.write_operator(
        filename::AbstractString,
        pack::WannierIO.SparseOperatorPack,
        ::WannierIO.JLD2Format;
        compress::Bool = _DEFAULT_JLD2_COMPRESS,
    )
    return JLD2.jldsave(filename, compress; pack)
end

function WannierIO.read_operator(filename::AbstractString, ::WannierIO.JLD2Format)
    return JLD2.load(filename, "pack")
end

end
