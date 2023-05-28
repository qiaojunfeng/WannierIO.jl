module WannierIO
using DocStringExtensions

include("common/const.jl")
include("common/type.jl")
include("util/fortran.jl")
include("util/lattice.jl")
include("util/header.jl")
include("util/toml.jl")
include("util/parser.jl")

using FortranFiles: FortranFile, FString, trimstring, Record

include("w90/win.jl")
include("w90/wout.jl")
include("w90/nnkp.jl")
include("w90/amn.jl")
include("w90/mmn.jl")
include("w90/eig.jl")
include("w90/chk.jl")
include("w90/unk.jl")
include("w90/spn.jl")
include("w90/band.jl")
include("w90/tb.jl")
include("w90/wsvec.jl")
include("w90/hr.jl")
include("w90/r.jl")
include("w90/hh_r.jl")

# volumetric files
include("volume/xsf.jl")
include("volume/cube.jl")
include("volume/bxsf.jl")

include("qe/band.jl")
include("qe/xml.jl")

end
