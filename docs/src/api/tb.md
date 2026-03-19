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
Pages   = ["operator/Rvector.jl"]
```

## Tight-binding operators

High-level TB API and conversions.

```@autodocs
Modules = [WannierIO]
Pages   = ["operator/tb.jl"]
```

## Hamiltonian operators

High-level HR API and conversions.

```@autodocs
Modules = [WannierIO]
Pages   = ["operator/hr.jl"]
```

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
