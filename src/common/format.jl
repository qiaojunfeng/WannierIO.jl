"""
Abstract supertype for all explicit file format tags used by `WannierIO`.

These tags are lightweight dispatch objects. Pass them to functions like
[`read_eig`](@ref), [`write_amn`](@ref), [`read_win`](@ref), etc. when you want
to choose a format explicitly instead of relying on auto detection.
"""
abstract type AbstractFileFormat end

"""
Abstract supertype for Fortran-style matrix/data file formats.
"""
abstract type AbstractFortranFormat <: AbstractFileFormat end

"""
Abstract supertype for Wannier90 small input-file formats such as `win` and `nnkp`.
"""
abstract type AbstractW90InputFormat <: AbstractFileFormat end

"""
Fortran formatted text IO.
"""
struct FortranText <: AbstractFortranFormat end

"""
Fortran unformatted record-based IO.
"""
struct FortranBinary <: AbstractFortranFormat end

"""
Fortran unformatted IO with stream access.

For example, file written using these Fortran code:
```fortran
OPEN(UNIT=11, FILE="ustream.demo", STATUS="NEW", ACCESS="STREAM", FORM="UNFORMATTED")
```
"""
struct FortranBinaryStream <: AbstractFortranFormat end

"""
Plain-text format for Wannier90 `win` and `nnkp` files.

The default `win` and `nnkp` files are text files, but they are not simple
tables that can be parsed with `readdlm`. This format tag selects the custom
Wannier90 text parser used by functions such as [`read_win`](@ref) and
[`read_nnkp`](@ref).

Compared with [`W90InputToml`](@ref), this format is for the native Wannier90
input syntax. Compared with [`FortranText`](@ref), it is for small structured
input files rather than large numeric matrix dumps such as `amn`, `mmn`, or `eig`.
"""
struct W90InputText <: AbstractW90InputFormat end

"""
TOML format for Wannier90 `win` and `nnkp` files.

This allows `win` or `nnkp` data to be stored in a standard TOML representation
that can be parsed by generic TOML tooling.
"""
struct W90InputToml <: AbstractW90InputFormat end

"""
    format_name(format)

Return a short human-readable name for a `FileFormat` tag.
"""
format_name(::FortranText) = "fortran-text"
format_name(::FortranBinary) = "fortran-binary"
format_name(::FortranBinaryStream) = "fortran-binary-stream"
format_name(::W90InputText) = "wannier90-text"
format_name(::W90InputToml) = "wannier90-toml"

Base.show(io::IO, format::AbstractFileFormat) = print(io, format_name(format))

"""
    fortran_format(; binary=false, stream=false)

Construct a Fortran-format tag.

- `binary=false` returns [`FortranText`](@ref)
- `binary=true, stream=false` returns [`FortranBinary`](@ref)
- `binary=true, stream=true` returns [`FortranBinaryStream`](@ref)
"""
function fortran_format(; binary::Bool=false, stream::Bool=false)
    binary || return FortranText()
    return stream ? FortranBinaryStream() : FortranBinary()
end

"""
    detect_fortran_format(file; stream=false)

Infer a Fortran-format tag from `file` using [`isbinary`](@ref).
"""
function detect_fortran_format(file::Union{IO,AbstractString}; stream::Bool=false)
    return isbinary(file) ? fortran_format(; binary=true, stream) : FortranText()
end

"""
    w90input_format(; toml=false)

Construct a format tag for `win`/`nnkp`-style input files.
"""
w90input_format(; toml::Bool=false) = toml ? W90InputToml() : W90InputText()

"""
    detect_w90input_format(file)

Infer a `win`/`nnkp` format tag from `file` using [`istoml`](@ref).
"""
function detect_w90input_format(file::Union{IO,AbstractString})
    return istoml(file) ? W90InputToml() : W90InputText()
end

"""
Native Wannier90 tight-binding text format (`.dat` files written by Wannier90).
"""
struct W90Dat <: AbstractFileFormat end

"""
HDF5 storage format for tight-binding data.

Requires loading the `HDF5` package.
"""
struct HDF5Format <: AbstractFileFormat end

"""
JLD2 storage format for tight-binding data.

Requires loading the `JLD2` package.
"""
struct JLD2Format <: AbstractFileFormat end

"""
Zarr storage format for tight-binding data.

Requires loading the `Zarr` package.
"""
struct ZarrFormat <: AbstractFileFormat end

format_name(::W90Dat) = "wannier90-dat"
format_name(::HDF5Format) = "hdf5"
format_name(::JLD2Format) = "jld2"
format_name(::ZarrFormat) = "zarr"

"""
    detect_w90dat_format(path)

Infer a tight-binding format tag from the file extension of `path`:
- `.h5` / `.hdf5` → [`HDF5Format`](@ref)
- `.jld2`          → [`JLD2Format`](@ref)
- `.zarr`          → [`ZarrFormat`](@ref)
    - anything else   → [`W90Dat`](@ref)
"""
function detect_w90dat_format(path::AbstractString)
    p = lowercase(path)
    if endswith(p, ".h5") || endswith(p, ".hdf5")
        return HDF5Format()
    elseif endswith(p, ".jld2")
        return JLD2Format()
    elseif endswith(p, ".zarr")
        return ZarrFormat()
    else
        return W90Dat()
    end
end
