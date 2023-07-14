# WannierIO.jl

A Julia package for reading/writing [Wannier90](https://github.com/wannier-developers/wannier90)
file formats.

The package is designed to be minimalistic to allow easy reuse in other packages.

## Quick examples

```julia
using WannierIO

A = read_amn("silicon.amn")
write_amn("silicon_2.amn", A)

chk = read_chk("silicon.chk")
write_chk("silicon_2.chk", chk; binary=true)
```

## Related packages

- [Wannier.jl](https://github.com/qiaojunfeng/Wannier.jl): Wannierization and Wannier interpolation.
    The IO part of `Wannier.jl` was isolated and moved into this package.
