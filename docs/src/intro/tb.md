# Tight-binding operators

```@meta
CurrentModule = WannierIO
```

This page describes the tight-binding (TB) files and the operator container
of WannierIO.

At a high level, WannierIO provides:

- reading and writing the Wannier90 `*_tb.dat`, `*_hr.dat`, `*_r.dat` and
  `*_wsvec.dat` files, exactly as stored,
- [`OperatorPack`](@ref), a container of named operators on a common set of
  R-vectors,
- sparse/compressed storage of an [`OperatorPack`](@ref) in HDF5, JLD2, Zarr
  and zipped Zarr files, see [Sparse storage](./sparse.md).

## Wannier90 files

Wannier90 TB files store each operator on the Wigner--Seitz R-vectors with
degeneracy weights (`Rdegens`), and `*_wsvec.dat` optionally adds the
translation vectors and degeneracies (`Tvectors`, `Tdegens`) of the
minimum-distance replica selection (MDRS).
WannierIO reads and writes these fields as they are, without applying them:

- [`TbDat`](@ref): `*_tb.dat`, the Hamiltonian and position operators,
- [`HrDat`](@ref): `*_hr.dat`, the Hamiltonian only (no lattice),
- [`RDat`](@ref): `*_r.dat`, the position operator only,
- [`WsvecDat`](@ref): `*_wsvec.dat`, with optional MDRS information.

```julia
using WannierIO

tbdat = read_w90_tb_dat("silicon_tb.dat")
wsvec = read_w90_wsvec_dat("silicon_wsvec.dat")

write_w90_tb_dat("silicon_copy_tb.dat", tbdat)
write_w90_wsvec_dat("silicon_copy_wsvec.dat", wsvec)
```

Absorbing the degeneracies and MDRS replicas into the operators is Wannier
interpolation, not file IO, and is done by
[Wannier.jl](https://github.com/qiaojunfeng/Wannier.jl): its
`read_w90_tb(prefix)` and `TightBindingModel(tbdat; wsvec)` return a model
whose operators are a plain inverse Fourier transform away from any kpoint.

## Operator container

An [`OperatorPack`](@ref) holds a header, the lattice, unique R-vectors and an
ordered dictionary of named operators. Each operator is an array
`n_wann × n_wann × … × n_Rvecs`: a scalar operator such as the Hamiltonian is
`n_wann × n_wann × n_Rvecs`, a vector operator such as the Berry connection is
`n_wann × n_wann × 3 × n_Rvecs`. The R-vectors carry no degeneracy weights, so
an operator at a kpoint is the plain sum
``O(\mathbf{k}) = \sum_{\mathbf{R}} e^{i 2\pi \mathbf{k} \cdot \mathbf{R}} O(\mathbf{R})``.

The container stores data and checks shapes; it does not interpret the
operators. Build one directly,

```julia
using WannierIO
using OrderedCollections: OrderedDict

operators = OrderedDict("hamiltonian" => H)  # n_wann × n_wann × n_Rvecs
pack = OperatorPack("my model", lattice, Rvectors, operators)
```

or export a Wannier.jl tight-binding model with `OperatorPack(tb)`, which
names the operators after the model's (`"hamiltonian"`,
`"berry_connection"`, …).

For all related type/function signatures, see the
[tight-binding API reference](../api/tb.md).
