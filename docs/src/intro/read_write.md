# Read and write files

```@meta
CurrentModule = WannierIO
```

WannierIO follows a consistent IO pattern across Wannier90 formats:

1. read with a `read_*` function,
2. if needed, modify the returned Julia data,
3. write back with the matching `write_*` function.

## IO conventions

- Function pairs are usually named `read_xxx` and `write_xxx`.
- Most readers/writers accept a filename or an `IO` stream.
- For several Wannier90 files, writing supports `binary=true` to emit Fortran unformatted files.
- Tight-binding files additionally support backend formats (HDF5/JLD2/Zarr) via
    [`read_w90_tb`](@ref) and [`write_w90_tb`](@ref).

## Common workflows

### 1. Round-trip a plain-text Wannier90 file

```julia
using WannierIO

amn = read_amn("silicon.amn")
write_amn("silicon_copy.amn", amn.A; header=amn.header)
```

The same pattern applies to other files such as
[`read_mmn`](@ref)/[`write_mmn`](@ref),
[`read_eig`](@ref)/[`write_eig`](@ref), and
[`read_nnkp`](@ref)/[`write_nnkp`](@ref).

### 2. Write Fortran binary files

```julia
using WannierIO

chk = read_chk("silicon.chk")
write_chk("silicon_binary.chk", chk; binary=true)

spn = read_spn("silicon.spn")
write_spn("silicon_binary.spn", spn; binary=true)
```

This is useful when interoperating with workflows that expect native Fortran
unformatted output.

### 3. Work with `win` input files

```julia
using WannierIO

win = read_win("silicon.win")
win["num_iter"] = 200
write_win("silicon_new.win", win)
```

`read_win` can parse both text and TOML-style inputs via format detection.

### 4. Read and write tight-binding files

```julia
using WannierIO

# High-level API with automatic wsvec handling when available.
pack = read_w90_tb("silicon_tb.dat")
write_w90_tb("silicon_out_tb.dat", pack)
```

For details on reduction and sparse backends, see
[Tight-binding operators](./tb.md) and [Sparse storage](./sparse.md).

### 5. Volumetric data formats

```julia
using WannierIO

xsf = read_xsf("charge.xsf")
write_xsf("charge_copy.xsf", xsf)

cube = read_cube("density.cube")
write_cube("density_copy.cube", cube)
```

Supported volumetric formats include XSF, CUBE, and BXSF.

## Tips

- Preserve metadata fields (for example headers/lattice metadata) when writing
    modified data.
- If you rely on binary files, test round-trip compatibility on your target
    compiler/platform because unformatted Fortran records are not universally portable.
- Prefer high-level helpers like [`read_w90_tb`](@ref) unless you explicitly need
    low-level file-to-file control.

For a deeper API-level reference, see [Wannier90](../api/w90.md) and
[Tight binding](../api/tb.md).
