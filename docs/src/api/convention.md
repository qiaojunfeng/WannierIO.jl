# Conventions

## Units

In most cases, the units of the function arguments and returned values are in
angstrom for lattice, and fractional (w.r.t. lattice) for atomic positions, etc.

## Variables

Here are some variable conventions used throughout the code and documentation.

### Names

- `U`: unitary transformation matrices
  - `A`: to differentiate between `U`, some times we use `A` specifically for
    the initial projection matrix, i.e., `amn` of Wannier90
- `M`: overlap matrices between neighboring kpoints, i.e., `mmn` of Wannier90

### Dimensions

Prefixed by `n_`, indicating this is an integer specifying the length of
some quantities; then followed by a short acronym of the quantity, to avoid
typing long names repeatedly.

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

## Functions

In most cases, there are some multiple dispatches for the same function.
For instance, there are three functions for reading the `chk` file:

```julia
read_chk(filename::AbstractString)
read_chk(filename::AbstractString, ::FortranText)
read_chk(filename::AbstractString, ::FortranBinary)
```

Why do we need several functions for reading the same format? The reasons are:

- Wannier90 is a Fortran code that accepts both Fortran text format and binary
  format (or called Fortran formatted and unformatted file, respectively).
  To avoid introducing lots of different function names for reading text or
  binary files, we use multiple dispatch to distinguish them: the type of the
  2nd argument (`FortranText` or `FortranBinary`) is used to distinguish the
  file type we want to read.
- However, asking the user to manually specify the file format is tedious,
  therefore, we provide a high-level user-friendly function (the 1st one), which
  - automatically detects the file format and calls the corresponding low-level
    (2nd or 3rd) function
  - prints some information of key quantities of the file, e.g., number of
    kpoints, number of Wannier functions, etc., to give user a hint of what
    have been parsed
  - hides some irrelevant return values, e.g., the header (the 1st line) of the
    file, since it has been printed in the stdout
- The low-level functions parse everything in the file, while the high-level
  function aims at user ergonomics.

Thus,

- In most cases using the high-level function is enough, e.g.,

  ```julia-repl
  julia> using WannierIO
  julia> A = read_amn("si2.amn")
  ┌ Info: Reading amn file
  │   filename = "si2.amn"
  │   header = "Created on  9Sep2022 at 16:41: 5"
  │   n_kpts = 8
  │   n_bands = 4
  └   n_wann = 4
  8-element Vector{Matrix{ComplexF64}}:
   [...]
  ```

- Use the low-level functions if you
  - want to get the header of the file (or quantities not returned by high-level function)
  - do not want stdout to be "polluted"

  ```julia-repl
  julia> using WannierIO
  julia> A, header = read_amn("si2.amn", WannierIO.FortranText())
  julia> typeof(A)
  Vector{Matrix{ComplexF64}} (alias for Array{Array{Complex{Float64}, 2}, 1})
  julia> header
  "Created on  9Sep2022 at 16:41: 5"
  ```

Note that usually the high-level function directly returns the quantities, e.g.,
a single `A` to avoid the user unpacking return values; however, often the low-level
functions return a `NamedTuple` of all the quantities, for the sake of clarity.

```julia-repl
julia> amn = read_amn("si2.amn", WannierIO.FortranText())
julia> typeof(amn)
NamedTuple{(:A, :header), Tuple{Vector{Matrix{ComplexF64}}, SubString{String}}}
```

Of course, when using low-level functions you can also directly access the
quantity without unpacking by

```julia-repl
julia> A = read_amn("si2.amn", WannierIO.FortranText()).A
```

When writing files, the user can specify whether to write in text or binary by
a keyword argument `binary` of the high-level function

```julia-repl
julia> write_amn("si2.amn", A; binary=true)
```

The `binary` keyword argument avoids the user specifying the file type
when calling the low-level functions, e.g.,

```julia-repl
julia> write_amn("si2.amn", A, WannierIO.FortranBinaryStream())
```
