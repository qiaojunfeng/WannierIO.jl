# Tight-binding models

```@meta
CurrentModule = WannierIO
```

API reference for tight-binding reduction and storage helpers.

For a narrative introduction and examples, see
[Tight-binding operators](../intro/tb.md) and [Sparse storage](../intro/sparse.md).

## R-vector reducers

Reducers that absorb R-degeneracy (and optional MDRS corrections) into operator
definitions.

```@autodocs
Modules = [WannierIO]
Pages   = ["util/Rvector.jl"]
```

## Operator containers

Container types and conversions (`TbDat`/`WsvecDat` <-> `OperatorPack`).

```@autodocs
Modules = [WannierIO]
Pages   = ["util/operator.jl"]
```

## Sparse representation

Sparse containers and conversion routines (`sparsify`, `densify`).

```@autodocs
Modules = [WannierIO]
Pages   = ["util/sparse.jl"]
```

## Page index

```@index
Modules = [WannierIO]
Pages   = ["tb.md"]
```
