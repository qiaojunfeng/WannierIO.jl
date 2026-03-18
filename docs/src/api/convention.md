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

Usually, the returned quantities are `Vector`s of some types (`Matrix{Float64}`,
`Vector{Float64}`, etc), and the indices follow the order of

1. kpoints
2. b-vectors (if needed)
3. bands
4. Wannier functions

For instance, the energy eigenvalues `eigenvalues` is a length-`n_kpts` vector,
with each element a length-`n_bands` vector of floats, i.e., `eigenvalues[ik][ib]`
is the `ib`-th eigenvalue at `ik`-th kpoint.

Here are some examples of indexing the vectors:

- `eigenvalues[ik][m]` for energy eigenvalues ``\varepsilon_{m \mathbf{k}}``
- `U[ik][m, n]` for the gauge matrix ``U_{mn \mathbf{k}}``
- `M[ik][ib][m, n]` for the overlap matrix ``M_{mn \mathbf{k}, \mathbf{k} + \mathbf{b}}``

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
