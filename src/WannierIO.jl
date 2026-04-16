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
include("w90/wsvec_dat.jl")
include("w90/tb_dat.jl")
include("w90/r_dat.jl")
include("w90/hh_r_dat.jl")
include("w90/u_mat.jl")
include("w90/isym.jl")

# volumetric files
include("volume/xsf.jl")
include("volume/cube.jl")
include("volume/bxsf.jl")

include("misc/epw.jl")

# Operators
include("operator/Rvector.jl")
include("operator/pack.jl")
include("operator/sparse.jl")
include("operator/tb.jl")
include("operator/hr.jl")

# include("precompile.jl")

end
