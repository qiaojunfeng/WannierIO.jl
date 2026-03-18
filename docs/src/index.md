# WannierIO.jl

A Julia package for reading/writing [Wannier90](https://github.com/wannier-developers/wannier90) files.

The package is designed to be minimalistic to allow easy reuse in other packages.

This package originates from the IO part of the
[Wannier.jl](https://github.com/qiaojunfeng/Wannier.jl) package.

## Wannier90 files

Input files: `amn`, `mmn`, `eig`, `chk`, `UNK`, `spn`, ...

Also support parsing both plain text and binary formats (in Fortran language,
called `formatted` and `unformatted` IO, respectively) for some files, e.g.,
`chk` and `UNK`.

Output files: `*_band.dat`, `*_tb.dat`, `*_wsvec.dat`, `*_hr.dat`, `*_r.dat`, `xsf`, `cube`, ...
