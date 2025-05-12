export read_uIu, write_uIu

"""
    read_uIu(filename)
    read_uIu(filename, ::FortranText; transpose_band_indices=true)
    read_uIu(filename, ::FortranBinary; transpose_band_indices=true)

Read the wannier90 `uIu` file.

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

function read_uIu(filename::AbstractString, ::FortranText; transpose_band_indices=true)
    return read_uHu(filename, FortranText(); transpose_band_indices)
end

function read_uIu(filename::AbstractString, ::FortranBinary; transpose_band_indices=true)
    return read_uHu(filename, FortranBinary(); transpose_band_indices)
end

function read_uIu(filename::AbstractString; kwargs...)
    if isbinary(filename)
        format = FortranBinary()
    else
        format = FortranText()
    end
    uIu, header = read_uIu(filename, format; kwargs...)

    n_kpts = length(uIu)
    @assert n_kpts > 0 "empty uIu matrix"
    n_bvecs = size(uIu[1], 1)
    @assert n_bvecs > 0 "empty uIu matrix"
    n_bands = size(uIu[1][1, 1], 1)
    @info "Reading uIu file" filename header n_kpts n_bvecs n_bands
    return uIu
end

"""
    write_uIu(filename, uIu; binary=false, header)
    write_uIu(filename, uIu, ::FortranText; header)
    write_uIu(filename, uIu, ::FortranBinary; header)

Write the `uIu` file.

# Keyword Arguments
- `transpose_band_indices`: see [`read_uIu`](@ref)
"""
function write_uIu end

function write_uIu(
    filename::AbstractString,
    uIu::AbstractVector,
    ::FortranText;
    header=default_header(),
    transpose_band_indices=true,
)
    return write_uHu(filename, uIu, FortranText(); header, transpose_band_indices)
end

function write_uIu(
    filename::AbstractString,
    uIu::AbstractVector,
    ::FortranBinary;
    header=default_header(),
    transpose_band_indices=true,
)
    return write_uHu(filename, uIu, FortranBinary(); header, transpose_band_indices)
end

function write_uIu(
    filename::AbstractString,
    uIu::AbstractVector;
    binary=false,
    header=default_header(),
    kwargs...,
)
    if binary
        format = FortranBinary()
    else
        format = FortranText()
    end
    n_kpts = length(uIu)
    @assert n_kpts > 0 "empty uIu matrix"
    n_bvecs = size(uIu[1], 1)
    @assert n_bvecs > 0 "empty uIu matrix"
    n_bands = size(uIu[1][1, 1], 1)
    @info "Writing uIu file" filename header n_kpts n_bands n_bvecs

    return write_uIu(filename, uIu, format; header, kwargs...)
end
