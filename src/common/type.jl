abstract type FileFormat end

"""
Fortran formatted IO.
"""
struct FortranText <: FileFormat end

"""
Fortran unformatted IO.
"""
struct FortranBinary <: FileFormat end

"""
Fortran unformatted IO with stream access.

For example, file written using these Fortran code:
```fortran
OPEN(UNIT=11, FILE="ustream.demo", STATUS="NEW", ACCESS="STREAM", FORM="UNFORMATTED")
```
"""
struct FortranBinaryStream <: FileFormat end

"""
Plain text format for Wannier90 `win` and `nnkp` files.

The W90 default `win` or `nnkp` are plain text files but are not
simple arrays of numbers that can be read by `readdlm`, therefore this struct
is used to indicate that the file is plain text but need to be handled
by corresponding functions, e.g., [`read_win`](@ref), [`read_nnkp`](@ref), etc.

This somewhat overlaps with [`FortranText`](@ref), but this one is only
used for small input parameter files e.g. `win` and `nnkp` (in comparison with
the [`Wannier90Toml`](@ref) format), while the [`FortranText`](@ref) (in
comparison with the [`FortranBinary`](@ref) format) is used for large matrices
e.g. `amn`, `mmn`, `eig`, etc.
"""
struct Wannier90Text <: FileFormat end

"""
TOML file format for Wannier90 `win` and `nnkp` files.

Here we introduce a TOML format for `win` and `nnkp`, so that once the `win` or
`nnkp` files are converted into TOML, the TOML files can be loaded by standard
TOML parsers without the headache of writing custom parsers in other Julia packages.

See also [`write_win`](@ref), [`write_nnkp`](@ref), etc.
"""
struct Wannier90Toml <: FileFormat end
