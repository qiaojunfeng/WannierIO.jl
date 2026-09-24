# Conventions

## Units

Unless explicitly stated otherwise:

- lattice vectors are in angstrom,
- reciprocal vectors are in angstrom^-1,
- fractional coordinates are expressed with respect to lattice vectors.

## Variables

The following naming conventions are used consistently across code and docs.

### Names

- `U`: unitary transformation matrices
  - `A`: initial projection matrix (`amn` in Wannier90 terminology)
- `M`: overlap matrices between neighboring kpoints, i.e., `mmn` of Wannier90

### Dimensions

Dimension variables are prefixed with `n_`.

- `n_bands`: number of bands
- `n_wann`: number of WFs
- `n_kpts`: number of kpoints
- `n_bvecs`: number of b-vectors
- `n_atoms`: number of atoms

### Indices

Band-resolved quantities are dense arrays with the kpoint as the last axis,
and the other axes in the order of

1. bands
2. Wannier functions
3. b-vectors (if needed)
4. kpoints

For instance, the energy eigenvalues `eigenvalues` is an `n_bands × n_kpts`
matrix, i.e., `eigenvalues[m, ik]` is the `m`-th eigenvalue at the `ik`-th
kpoint.

Here are some examples of indexing the arrays:

- `eigenvalues[m, ik]` for energy eigenvalues ``\varepsilon_{m \mathbf{k}}``
- `U[m, n, ik]` for the gauge matrix ``U_{mn \mathbf{k}}``
- `M[m, n, ib, ik]` for the overlap matrix ``M_{mn \mathbf{k}, \mathbf{k} + \mathbf{b}}``

where

- `ik`: index of kpoints
- `ib`: index of b-vectors
- `m`: index of bands
- `n`: index of Wannier functions

### Containers

Reader/writer APIs follow a simple rule for returned/accepted grouped data:

- If a parser returns up to 3 values, it returns a `NamedTuple`
- If a parser returns more than 3 values, it returns a thin container `struct`

The same thin container structs are accepted by corresponding writer functions.
This keeps small APIs lightweight while giving large file formats a centralized
data model that downstream packages can reuse.

## Functions

Most top-level APIs have multiple dispatch variants.
For example, there are format-specific methods for reading `chk` files:

```julia
read_chk(filename::AbstractString)
read_chk(filename::AbstractString, ::FortranText)
read_chk(filename::AbstractString, ::FortranBinary)
```

Why this design:

- high-level functions are convenient and automatically detect formats,
- low-level methods remain available when explicit format control is needed.

Thus,

- In most cases, use the high-level function:

  ```julia-repl
  julia> using WannierIO
  julia> amn = read_amn("si2.amn");
  julia> amn.header
  "Created on  9Sep2022 at 16:41: 5"
  julia> amn.A
  8-element Vector{Matrix{ComplexF64}}:
   [...]
  ```

- Use low-level methods when you need explicit format control.

  ```julia-repl
  julia> using WannierIO
  julia> amn = read_amn("si2.amn", WannierIO.FortranText())
  julia> amn.header
  "Created on  9Sep2022 at 16:41: 5"
  ```

When writing files, high-level methods commonly expose a `binary` keyword:

```julia-repl
julia> write_amn("si2.amn", A; binary=true)
```

This avoids the need to call format-specific low-level methods in typical use.
