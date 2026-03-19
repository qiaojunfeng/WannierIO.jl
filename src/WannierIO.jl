"""
`WannierIO.jl`: a package for reading and writing Wannier90 file formats.

---

$(README)

---

Exported functions:

$(EXPORTS)

"""
module WannierIO

using Printf: @printf, @sprintf
using DocStringExtensions

using LinearAlgebra
using SparseArrays
using StaticArrays

using CrystalBase

include("common/const.jl")
include("common/format.jl")
include("util/fortran.jl")
include("util/header.jl")
include("util/toml.jl")
include("util/parser.jl")
include("util/compare.jl")
include("util/Rvector.jl")
include("util/operator.jl")
include("util/sparse.jl")

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
include("w90/uHu.jl")
include("w90/uIu.jl")
include("w90/band.jl")
include("w90/hr_dat.jl")
include("w90/hr.jl")
include("w90/wsvec.jl")
include("w90/tb_dat.jl")
include("w90/tb.jl")
include("w90/r.jl")
include("w90/hh_r.jl")
include("w90/u_mat.jl")
include("w90/isym.jl")

# volumetric files
include("volume/xsf.jl")
include("volume/cube.jl")
include("volume/bxsf.jl")

include("misc/epw.jl")

end
