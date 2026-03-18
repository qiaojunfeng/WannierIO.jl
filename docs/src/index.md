# WannierIO.jl

A Julia package for reading and writing
[Wannier90](https://github.com/wannier-developers/wannier90) file formats.

The package is intentionally minimal so the IO layer can be reused in other
Wannier and electronic-structure tooling.

This package originates from the IO part of the
[Wannier.jl](https://github.com/qiaojunfeng/Wannier.jl) package.

## What you can do

- Read/write core Wannier90 formats (`win`, `amn`, `mmn`, `eig`, `chk`, `UNK`, `spn`, ...)
- Read/write tight-binding datasets (`*_tb.dat`, `*_wsvec.dat`, `*_hr.dat`, `*_r.dat`)
- Read/write volumetric formats (`xsf`, `cube`, `bxsf`)
- Store tight-binding operators in sparse/compressed backends (HDF5/JLD2/Zarr)

## Quick start

```julia
using WannierIO

amn = read_amn("silicon.amn")
write_amn("silicon_copy.amn", amn.A; header=amn.header)

tb = read_w90_tb("silicon_tb.dat")
write_w90_tb("silicon_out_tb.dat", tb)
```

## Documentation map

- Introduction: read/write workflows, tight-binding reduction, sparse storage.
- API reference: conventions, utility types, Wannier90/TB APIs, volumetric and EPW helpers.
