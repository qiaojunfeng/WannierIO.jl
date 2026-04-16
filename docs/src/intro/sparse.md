# Sparse storage for tight-binding operators

```@meta
CurrentModule = WannierIO
```

Wannier90 tight-binding datasets can become large because each R-vector stores
full dense operator matrices (`H`, `rx`, `ry`, `rz`).

WannierIO supports a sparse-packed representation that:

- drops near-zero matrix elements using a tolerance,
- stores nonzeros in compact CSC-based packs,
- writes compressed binary formats for efficient disk storage.

## Supported backends

Sparse TB IO is exposed through [`read_operator`](@ref) and [`write_operator`](@ref)
with format objects:

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

# Read Wannier90 tb.dat (+ wsvec when available).
pack = read_w90_tb("Si2_valence_tb.dat")

# Write sparse/compressed HDF5.
write_operator("Si2_valence.h5", pack, WannierIO.HDF5Format(); atol=1e-10)

# Read back (returns dense OperatorPack).
pack_h5 = read_operator("Si2_valence.h5")
```

If no format is provided, [`detect_operator_format`](@ref) selects it from the
filename extension when calling [`read_operator`](@ref) or [`write_operator`](@ref).

## Precision and sparsification controls

Sparsification behavior is controlled by [`SparseOption`](@ref), usually passed
as keyword arguments to `write_operator`:

- `atol`: threshold below which real/imaginary parts are dropped,
- `value_type`: numeric type for stored values (for example `Float32`),
- `index_type`: integer type for sparse indices (for example `Int16`).

```julia
using WannierIO
using JLD2

pack = read_w90_tb("Si2_valence_tb.dat")
write_operator(
    "Si2_valence.jld2",
    pack,
    WannierIO.JLD2Format();
    atol=1e-10,
    value_type=Float32,
    index_type=Int16,
)
```

## Backend-specific compression options

In addition to `atol/value_type/index_type`, each backend accepts extra controls:

- `HDF5Format`: `deflate`, `shuffle`
- `JLD2Format`: `compress`
- `ZarrFormat`/`ZarrZipFormat`: `clevel`, `shuffle`

```julia
using WannierIO
using Zarr

pack = read_w90_tb("Si2_valence_tb.dat")
write_operator("Si2_valence.zarr", pack, ZarrFormat(); clevel=7, shuffle=true)
```

## Direct sparse conversion API

For workflows that do not immediately write files, you can convert explicitly:

```julia
using WannierIO

pack = read_w90_tb("Si2_valence_tb.dat")
spack = sparsify(pack; atol=1e-10, value_type=Float32)
pack2 = densify(spack)
```

The sparse data container is [`SparseOperatorPack`](@ref), which stores each
operator as a [`CscPack`](@ref).

## Notes

- Reading sparse backends returns dense [`OperatorPack`](@ref) for easy use in
    interpolation and downstream code.
- Choose `atol` carefully: larger values reduce file size more aggressively,
    but can remove physically relevant small matrix elements.
