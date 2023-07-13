# Conventions

## Units

In most cases, the units of the function arguments and returned values are in
angstrom for lattice, and fractional (w.r.t. lattice) for atomic positions, etc.

## Variable names

Here we list the common variable names used throughout the code and
documentation. If you want to introduce a new variable that might be used
widely, consider following the same style.

### length of quantities

Prefixed by `n_`, indicating this is an integer specifying the length of
some quantities; then followed by a short acronym of the quantity, to avoid
typing long names repeatedly.

- `n_bands`: number of bands
- `n_wann`: number of WFs
- `n_kpts`: number of kpoints
- `n_bvecs`: number of b-vectors
- `n_atoms`: number of atoms

### Matrices

- `U`: unitary transformation matrices
  - `A`: to differentiate between `U`, some times we use `A` specifically for
    the initial projection matrix, i.e., `amn` of Wannier90
- `M`: overlap matrices between neighboring kpoints, i.e., `mmn` of Wannier90
- `E`: energy eigenvalues, i.e., `eig` of Wannier90
- `S`: spin operator, i.e., `spn` of Wannier90
- `R`: the R-vectors used in Wannier interpolation

## Indices

Usually, the returned quantities are `Vector`s of some types (`Matrix{Float64}`,
`Vector{Float64}`, etc), and the indices follow the order of

1. kpoints
2. b-vectors (if needed)
3. bands
4. Wannier functions

For instance, the energy eigenvalues `E` is a length-`n_kpts` vector, with each
element a length-`n_bands` vector of floats, i.e., `E[ik][ib]` is the `ib`-th
eigenvalue at `ik`-th kpoint.

Here are some examples of indexing the vectors:

- `E[ik][m]`
- `U[ik][m, n]`
- `M[ik][ib][m, n]`

where

- `ik`: index of kpoints
- `ib`: index of b-vectors
- `m`: index of bands
- `n`: index of Wannier functions
