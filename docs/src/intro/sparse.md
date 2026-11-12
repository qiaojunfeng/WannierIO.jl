# Sparse storage for tight-binding operators

```@meta
CurrentModule = WannierIO
```

Tight-binding operators can become large because each R-vector stores full
dense operator matrices, for the Hamiltonian and for every component of the
Berry connection or the spin.

A [`SparseOperatorPack`](@ref) stores the operators of an
[`OperatorPack`](@ref) as compressed sparse column matrices ([`CscPack`](@ref)),
and is what the binary formats write:

- convert with `SparseOperatorPack(pack)`, choosing the precision there,
- write and read the sparse pack with [`write_operator`](@ref) and
  [`read_operator`](@ref), which store it as it is,
- convert back with `OperatorPack(spack)`.

## Supported backends

[`write_operator`](@ref) and [`read_operator`](@ref) take a format object:

- [`HDF5Format`](@ref)
- [`JLD2Format`](@ref)
- [`ZarrFormat`](@ref)
- [`ZarrZipFormat`](@ref)

These methods are provided by optional extensions. Load the corresponding package
before using the format.

```julia
using Pkg
Pkg.add(["HDF5", "JLD2", "Zarr"])

using WannierIO
using HDF5   # or JLD2 / Zarr
```

## Quick start

```julia
using WannierIO
using HDF5

# `pack` is an OperatorPack, e.g. `OperatorPack(tb)` of a Wannier.jl model.
write_operator("Si2_valence.h5", SparseOperatorPack(pack))

# Read back the SparseOperatorPack, and its dense operators.
pack_h5 = OperatorPack(read_operator("Si2_valence.h5"))
```

If no format is provided, [`detect_operator_format`](@ref) selects it from the
filename extension when calling [`read_operator`](@ref) or [`write_operator`](@ref).

## Precision

The conversion to a [`SparseOperatorPack`](@ref) is lossless by default: it
keeps the precision of the pack and drops only the elements that are exactly
zero. Two choices make the file smaller, at a cost in accuracy:

- the value type `Tv`: `SparseOperatorPack{Float32}(pack)` rounds each element
  to about 7 significant digits,
- `atol`: drops the real and imaginary parts at most `atol` in magnitude. It is
  absolute, in the units of each operator (eV for a Hamiltonian, Å for a Berry
  connection).

An operator whose imaginary parts are all at most `atol` is stored as real.
For a 70-Wannier-function Cu₂P₂Zr model with spin-orbit coupling (923
R-vectors, Hamiltonian and Berry connection; `tb.dat` and `wsvec.dat` take
950 MB), the HDF5 file is

| conversion                                       | size   | max change in the bands |
|:-------------------------------------------------|:-------|:------------------------|
| `SparseOperatorPack(pack)`                       | 110 MB | 0                       |
| `SparseOperatorPack{Float32}(pack)`              | 72 MB  | 0.8 µeV                 |
| `SparseOperatorPack{Float32}(pack; atol = 1e-6)` | 45 MB  | 16 µeV                  |

The cost of `atol` depends on the model. Across the 20 `tb.dat` models of
WannierDatasets (4 to 26 Wannier functions), all three conversions read back
as written, the lossless file is about 3.5 times smaller than `tb.dat` and
`wsvec.dat`, and `Float32` halves it again with bands within about 1 µeV. The
further gain of `atol = 1e-6` ranges widely:

| model                     | lossless | `Float32`          | `Float32`, `atol = 1e-6` |
|:--------------------------|:---------|:-------------------|:-------------------------|
| Fe, collinear (9 WFs)     | 1.1 MB   | 0.6 MB, 0.7 µeV    | 0.6 MB, 7 µeV            |
| GaAs (8 WFs)              | 3.0 MB   | 1.5 MB, 0.4 µeV    | 1.3 MB, 34 µeV           |
| TiO₂ (26 WFs)             | 12.2 MB  | 6.1 MB, 0.3 µeV    | 3.1 MB, 23 µeV           |
| Si, `3d` (18 WFs)         | 14.8 MB  | 7.3 MB, 0.3 µeV    | 3.8 MB, 130 µeV          |

Rounding to `Float32` costs little, since each element keeps its relative
precision. `atol` costs more: it removes every part at most `atol` in
magnitude, and in Wannier90 models most of these are not whole elements but
the small imaginary parts of kept ones, whose sum over the R-vectors moves the
bands. `Float32` is therefore the safe compact choice; check the bands of the
model before relying on an `atol`. Larger values of `atol` reduce the file
further, but can remove physically relevant small matrix elements.

The conversion to real rarely applies to Wannier90 output: even maximally
localized functions that are real in principle carry imaginary parts up to
`1e-4` (Si) in their operators, far above the `atol` of a compact file. It helps for operators that are real by
construction.

The index type `Ti` (default `Int32`) also counts the nonzeros of each
operator, so a smaller type overflows for all but the smallest models.

## Backend-specific compression options

The keywords of [`write_operator`](@ref) are the compression options of each
backend:

- `HDF5Format`: `deflate`, `shuffle`
- `JLD2Format`: `compress`
- `ZarrFormat`/`ZarrZipFormat`: `clevel`, `shuffle`

```julia
using WannierIO
using Zarr

write_operator("Si2_valence.zarr", SparseOperatorPack(pack), ZarrFormat(); clevel=7, shuffle=true)
```

## Sparse containers

Each operator of a [`SparseOperatorPack`](@ref) is a [`CscPack`](@ref): its
`n_wann × n_wann` matrices, one per entry of the trailing (component and
R-vector) axes in column-major order, in compressed sparse column format,
together with the `dims` of the dense array. `CscPack(op; atol)` packs one
dense operator and `Array(csc)` unpacks it.
