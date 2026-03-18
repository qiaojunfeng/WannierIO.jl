# Utility

Shared utility types and helper functions used across file-format readers,
writers, and conversion pipelines.

```@meta
CurrentModule = WannierIO
```

## Constants

Physical constants and package-wide defaults.

```@autodocs
Modules = [WannierIO]
Pages   = ["common/const.jl"]
```

## Helper functions

General helpers for parsing, headers, TOML conversion, and comparisons.

```@autodocs
Modules = [WannierIO]
Pages   = ["util/header.jl", "util/toml.jl", "util/parser.jl", "util/compare.jl"]
```

## Fortran related

Low-level helpers for formatted/unformatted Fortran IO.

```@autodocs
Modules = [WannierIO]
Pages   = ["util/fortran.jl"]
```

## File formats

Format tags used by high-level APIs for dispatch and auto-detection.

```@autodocs
Modules = [WannierIO]
Pages   = ["common/format.jl"]
```

## Package docstring

Top-level package docstring and exported-symbol overview.

```@autodocs
Modules = [WannierIO]
Pages   = ["WannierIO.jl"]
```

## Page index

```@index
Modules = [WannierIO]
Pages   = ["util.md"]
```
