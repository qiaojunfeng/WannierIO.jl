export read_uIu, write_uIu

"""
    read_uIu(filename)
    read_uIu(file, ::FortranText; transpose_band_indices=true)
    read_uIu(file, ::FortranBinary; transpose_band_indices=true)

Read the wannier90 `uIu` file.

# Arguments
- `file`: The name of the input file, or an `IO`.

# Keyword Arguments
- `transpose_band_indices`: QE pw2wannier90.x writes the matrix in a strange
    transposed manner; if reading a QE-generated `uIu` file, this flag should
    be true to restore the band indices order, so that the returned matrix
    has the correct order, i.e.,
    `uIu[ik][ib1, ib2][m, n]` is
    ``\\langle u_{m, k + b_1} | u_{n, k + b_2} \\rangle``

# Return
- `uIu`: a length-`n_kpts` vector, each element is a `n_bvecs * n_bvecs` matrix,
    then each element is a `n_bands * n_bands` matrix
- `header`: 1st line of the file
"""
function read_uIu end

function read_uIu(filename::AbstractString, ::FortranText; transpose_band_indices = true)
    return read_uHu(filename, FortranText(); transpose_band_indices)
end

function read_uIu(io::IO, ::FortranText; transpose_band_indices = true)
    return read_uHu(io, FortranText(); transpose_band_indices)
end

function read_uIu(filename::AbstractString, ::FortranBinary; transpose_band_indices = true)
    return read_uHu(filename, FortranBinary(); transpose_band_indices)
end

function read_uIu(io::FortranFile, ::FortranBinary; transpose_band_indices = true)
    return read_uHu(io, FortranBinary(); transpose_band_indices)
end

function read_uIu(filename::AbstractString; kwargs...)
    format = detect_fortran_format(filename)
    return read_uIu(filename, format; kwargs...)
end

"""
    write_uIu(filename, uIu; binary=false, header)
    write_uIu(file, uIu, ::FortranText; header)
    write_uIu(file, uIu, ::FortranBinary; header)

Write the `uIu` file.

# Arguments
- `file`: The name of the output file, or an `IO`.
- `uIu`: a length-`n_kpts` vector, each element is a `n_bvecs * n_bvecs` matrix,
    then each element is a `n_bands * n_bands` matrix

# Keyword Arguments
- `transpose_band_indices`: see [`read_uIu`](@ref)
"""
function write_uIu end

function write_uIu(
        filename::AbstractString,
        uIu::AbstractVector,
        ::FortranText;
        header = default_header(),
        transpose_band_indices = true,
    )
    return write_uHu(filename, uIu, FortranText(); header, transpose_band_indices)
end

function write_uIu(
        io::IO,
        uIu::AbstractVector,
        ::FortranText;
        header = default_header(),
        transpose_band_indices = true,
    )
    return write_uHu(io, uIu, FortranText(); header, transpose_band_indices)
end

function write_uIu(
        filename::AbstractString,
        uIu::AbstractVector,
        ::FortranBinary;
        header = default_header(),
        transpose_band_indices = true,
    )
    return write_uHu(filename, uIu, FortranBinary(); header, transpose_band_indices)
end

function write_uIu(
        io::FortranFile,
        uIu::AbstractVector,
        ::FortranBinary;
        header = default_header(),
        transpose_band_indices = true,
    )
    return write_uHu(io, uIu, FortranBinary(); header, transpose_band_indices)
end

function write_uIu(
        filename::AbstractString,
        uIu::AbstractVector;
        binary = false,
        header = default_header(),
        kwargs...,
    )
    format = fortran_format(; binary)
    return write_uIu(filename, uIu, format; header, kwargs...)
end
