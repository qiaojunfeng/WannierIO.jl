export read_eig, write_eig

"""
    $(SIGNATURES)

Reshape a vector of eigenvalues into a matrix of eigenvalues.

Auto detect the number of bands and kpoints.
"""
@inline function _reshape_eig(
    idx_b::AbstractVector, idx_k::AbstractVector, eig::AbstractVector
)
    # find unique elements
    n_bands = length(Set(idx_b))
    n_kpts = length(Set(idx_k))
    E = reshape(eig, (n_bands, n_kpts))

    return map(ik -> E[:, ik], 1:n_kpts)
end

"""
    $(SIGNATURES)

Check that eigenvalues are in order.

Some times there are small noises, use `digits` to set the number of digits for comparisons.
"""
@inline function _check_eig_order(eigenvalues::AbstractVector; digits=7)
    round_digits(x) = round(x; digits)
    for (ik, eig) in enumerate(eigenvalues)
        if !issorted(eig; by=round_digits)
            @warn "Eigenvalues are not sorted at " ik eig
        end
    end
end

"""
    read_eig(file)
    read_eig(file, ::FortranText)
    read_eig(file, ::FortranBinaryStream)

Read the wannier90 `eig` file.

# Arguments
- `file`: The name of the input file, or an `IO`.

# Return
- `eigenvalues`: a lenth-`n_kpts` vector, each element is a length-`n_bands` vector

The 1st version is a convenience wrapper function for the 2nd and 3rd versions.
"""
function read_eig end

function read_eig(io::IO, ::FortranText)
    lines = readlines(io)

    n_lines = length(lines)
    idx_b = zeros(Int, n_lines)
    idx_k = zeros(Int, n_lines)
    eig = zeros(Float64, n_lines)

    for i in 1:n_lines
        arr = split(lines[i])
        idx_b[i] = parse(Int, arr[1])
        idx_k[i] = parse(Int, arr[2])
        eig[i] = parse(Float64, arr[3])
    end

    eigenvalues = _reshape_eig(idx_b, idx_k, eig)
    _check_eig_order(eigenvalues)
    return eigenvalues
end

function read_eig(io::IO, ::FortranBinaryStream)
    idx_b = Vector{Int}()
    idx_k = Vector{Int}()
    eig = Vector{Float64}()

    # gfortran integer is 4 bytes
    Tint = Int32

    while !eof(io)
        push!(idx_b, read(io, Tint))
        push!(idx_k, read(io, Tint))
        push!(eig, read(io, Float64))
    end

    eigenvalues = _reshape_eig(idx_b, idx_k, eig)
    _check_eig_order(eigenvalues)
    return eigenvalues
end

function read_eig(filename::AbstractString, format::FileFormat)
    return open(filename) do io
        read_eig(io, format)
    end
end

function read_eig(file::Union{IO,AbstractString})
    format = isbinary(file) ? FortranBinaryStream() : FortranText()
    return read_eig(file, format)
end

"""
    write_eig(file, eigenvalues; binary=false)
    write_eig(file, eigenvalues, ::FortranText)
    write_eig(file, eigenvalues, ::FortranBinaryStream)

Write `eig` file.

# Arguments
- `file`: The name of the output file, or an `IO`.
- `eigenvalues`: a length-`n_kpts` vector, each element is a length-`n_bands` vector

# Keyword arguments
- `binary`: if true write in Fortran binary format.
"""
function write_eig end

function write_eig(io::IO, eigenvalues::AbstractVector, ::FortranText)
    n_kpts = length(eigenvalues)
    n_kpts > 0 || throw(ArgumentError("Empty eigenvalues"))
    n_bands = length(eigenvalues[1])

    for ik in 1:n_kpts
        for ib in 1:n_bands
            @printf(io, "%5d%5d%18.12f\n", ib, ik, eigenvalues[ik][ib])
        end
    end
end

function write_eig(io::IO, eigenvalues::AbstractVector, ::FortranBinaryStream)
    n_kpts = length(eigenvalues)
    n_kpts > 0 || throw(ArgumentError("Empty eigenvalues"))
    n_bands = length(eigenvalues[1])

    # gfortran integer is 4 bytes
    Tint = Int32
    # I write in Fortran stream IO, so just plain julia `open`
    for ik in 1:n_kpts
        for ib in 1:n_bands
            write(io, Tint(ib))
            write(io, Tint(ik))
            write(io, Float64(eigenvalues[ik][ib]))
        end
    end
end

function write_eig(
    filename::AbstractString, eigenvalues::AbstractVector, format::FileFormat
)
    open(filename, "w") do io
        write_eig(io, eigenvalues, format)
    end
end

function write_eig(
    file::Union{IO,AbstractString}, eigenvalues::AbstractVector; binary=false
)
    format = binary ? FortranBinaryStream() : FortranText()
    write_eig(file, eigenvalues, format)
end
