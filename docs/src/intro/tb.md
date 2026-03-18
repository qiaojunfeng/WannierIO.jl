# Tight-binding models

```@meta
CurrentModule = WannierIO
```

This page describes the tight-binding (TB) workflow in WannierIO.

At a high level, the TB API provides:

- reading and writing native Wannier90 `*_tb.dat` and `*_wsvec.dat` files,
- automatic reduction of R-vector degeneracies and MDRS corrections,
- conversion to a clean operator container for interpolation and post-processing,
- optional sparse/compressed storage backends (HDF5/JLD2/Zarr).

## R-vector reduction notes

Wannier90 TB files may include R-degeneracies (the Wigner-Seitz method),
or, T-vectors and degeneracies (the minima distance replica selection (MDRS) method).
In Wannier90, these are stored in the `*_tb.dat` and `*_wsvec.dat` files.

WannierIO handles these through [`RvectorReducer`](@ref), so the resulting
[`OperatorPack`](@ref) is expressed on a reduced set of R vectors, and the
degeneracies are absorbed into the operator definitions.

This means downstream interpolation can be implemented as a straightforward
inverse Fourier transform on the reduced operators.

## Core data flow

The key data structures are:

- [`TbDat`](@ref): dense data from `*_tb.dat`, including Hamiltonian and position operators,
- [`WsvecDat`](@ref): data from `*_wsvec.dat` with optional MDRS information,
- [`OperatorPack`](@ref): reduced operator representation used by the high-level TB API,
- [`SparseOperatorPack`](@ref): sparse-packed representation for compact storage.

For most users, [`read_w90_tb`](@ref) is the recommended entry point.
It loads `*_tb.dat`, auto-detects and uses the matching `*_wsvec.dat` when available,
and returns a reduced [`OperatorPack`](@ref) that is ready for further processing.

## Typical usage

### 1. Read and write native Wannier90 TB files

```julia
using WannierIO

# Reads prefix_tb.dat, and if present, prefix_wsvec.dat.
pack = read_w90_tb("silicon_tb.dat")

# Write a reduced OperatorPack back to a native tb.dat file.
write_w90_tb("silicon_reduced_tb.dat", pack)
```

### 2. Keep explicit control of wsvec handling

```julia
using WannierIO

tb = read_w90_tb_dat("silicon_tb.dat")
ws = read_w90_wsvec("silicon_wsvec.dat")

# Apply reduction explicitly.
pack = pack(tb, ws)

# Write both files explicitly.
write_w90_tb("silicon_out_tb.dat", tb, ws)
```

### 3. Save in sparse/compressed formats

```julia
using WannierIO
using HDF5   # or JLD2 / Zarr

pack = read_w90_tb("silicon_tb.dat")
write_w90_tb("silicon.h5", pack, HDF5Format(); atol=1e-10)

pack2 = read_w90_tb("silicon.h5", HDF5Format())
```

The same high-level API works for [`JLD2Format`](@ref) and [`ZarrFormat`](@ref).
See [sparse storage](./sparse.md) for backend-specific details and precision options.
