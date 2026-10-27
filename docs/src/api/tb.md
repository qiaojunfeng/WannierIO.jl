# Tight-binding models

```@meta
CurrentModule = WannierIO
```

API reference for the tight-binding operator container and its storage.

For a narrative introduction and examples, see
[Tight-binding operators](../intro/tb.md) and [Sparse storage](../intro/sparse.md).

## Operator containers

Container types and conversions for generic operator storage.

```@autodocs
Modules = [WannierIO]
Pages   = ["operator/pack.jl"]
```

## Sparse representation

Sparse containers and conversion routines (`sparsify`, `densify`).

```@autodocs
Modules = [WannierIO]
Pages   = ["operator/sparse.jl"]
```

## Page index

```@index
Modules = [WannierIO]
Pages   = ["tb.md"]
```
