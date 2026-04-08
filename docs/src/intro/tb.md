# Tight-binding models

```@meta
CurrentModule = WannierIO
```

This page describes the tight-binding (TB) workflow in WannierIO.

At a high level, the TB API provides:

- reading and writing Wannier90 `*_tb.dat` and `*_wsvec.dat` files,
- automatic reduction of R-vector degeneracies and MDRS corrections,
- conversion to a clean operator container for interpolation and post-processing,
- optional sparse/compressed storage backends (HDF5/JLD2/Zarr/ZarrZip).

## R-vector reduction notes

Wannier90 TB files may include R-degeneracies (Wigner-Seitz method), and
optionally T-vectors/degeneracies from the minima distance replica selection
(MDRS) method.
In Wannier90, these are stored in the `*_tb.dat` and `*_wsvec.dat` files.

WannierIO handles these through [`RvectorReducer`](@ref), so the resulting
[`OperatorPack`](@ref) is expressed on a reduced set of R vectors, and the
degeneracies are absorbed into the operator definitions.

This means downstream interpolation can be implemented as a straightforward
inverse Fourier transform on the reduced operators.

## Core data flow

Key data structures:

- [`TbDat`](@ref): dense data from `*_tb.dat`, including Hamiltonian and position operators,
- [`WsvecDat`](@ref): data from `*_wsvec.dat` with optional MDRS information,
- [`OperatorPack`](@ref): reduced operator representation used by the high-level TB API,
- [`SparseOperatorPack`](@ref): sparse-packed representation for compact storage.

For most users, [`read_w90_tb`](@ref) is the recommended entry point.
It loads `*_tb.dat`, auto-detects and uses the matching `*_wsvec.dat` when available,
and returns a reduced [`OperatorPack`](@ref) that is ready for further processing.

This `read_w90_tb`/`write_w90_tb` API is for Wannier90 text files.
For sparse/compressed backends, use [`read_operator`](@ref) and [`write_operator`](@ref)
with [`HDF5Format`](@ref), [`JLD2Format`](@ref), [`ZarrFormat`](@ref),
or [`ZarrZipFormat`](@ref).

## Typical usage

### 1. Read and write Wannier90 TB files

```julia
using WannierIO

# Reads prefix_tb.dat, and if present, prefix_wsvec.dat.
pack = read_w90_tb("silicon_tb.dat")

# Write an OperatorPack back to a tb.dat file.
write_w90_tb("silicon_reduced_tb.dat", pack)
```

### 2. Keep explicit control of wsvec handling

```julia
using WannierIO

tb = read_w90_tb_dat("silicon_tb.dat")
ws = read_w90_wsvec_dat("silicon_wsvec.dat")

# Apply reduction explicitly.
pack = pack(tb, ws)

# Write both files explicitly.
write_w90_tb("silicon_out_tb.dat", tb, ws)
```

For all related type/function signatures, see the
[tight-binding API reference](../api/tb.md).
