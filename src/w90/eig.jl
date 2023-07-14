using Printf: @printf

export read_eig, write_eig

"""
Reshape a vector of eigenvalues into a matrix of eigenvalues.
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
Check that eigenvalues are in order, some times there are small noises.
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
    read_eig(filename)
    read_eig(filename, ::FortranText)
    read_eig(filename, ::FortranBinaryStream)

Read the wannier90 `eig` file.

# Return
- `eigenvalues`: a lenth-`n_kpts` vector, each element is a length-`n_bands` vector

The 1st version is a convenience wrapper function for the 2nd and 3rd versions.
"""
function read_eig end

function read_eig(filename::AbstractString, ::FortranText)
    lines = open(filename) do io
        readlines(io)
    end

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

function read_eig(filename::AbstractString, ::FortranBinaryStream)
    idx_b = Vector{Int}()
    idx_k = Vector{Int}()
    eig = Vector{Float64}()

    # gfortran integer is 4 bytes
    Tint = Int32

    open(filename) do io
        while !eof(io)
            push!(idx_b, read(io, Tint))
            push!(idx_k, read(io, Tint))
            push!(eig, read(io, Float64))
        end
    end

    eigenvalues = _reshape_eig(idx_b, idx_k, eig)
    _check_eig_order(eigenvalues)
    return eigenvalues
end

function read_eig(filename::AbstractString)
    if isbinary(filename)
        format = FortranBinaryStream()
    else
        format = FortranText()
    end
    eigenvalues = read_eig(filename, format)

    n_kpts = length(eigenvalues)
    @assert n_kpts > 0 "Empty eig file"
    n_bands = length(eigenvalues[1])
    @info "Reading eig file" filename n_kpts n_bands

    return eigenvalues
end

"""
    write_eig(filename, eigenvalues; binary=false)
    write_eig(filename, eigenvalues, ::FortranText)
    write_eig(filename, eigenvalues, ::FortranBinaryStream)

Write `eig` file.

# Arguments
- `eigenvalues`: a length-`n_kpts` vector, each element is a length-`n_bands` vector

# Keyword arguments
- `binary`: if true write in Fortran binary format.
"""
function write_eig end

function write_eig(filename::AbstractString, eigenvalues::AbstractVector, ::FortranText)
    n_kpts = length(eigenvalues)
    @assert n_kpts > 0 "Empty eigenvalues"
    n_bands = length(eigenvalues[1])

    open(filename, "w") do io
        for ik in 1:n_kpts
            for ib in 1:n_bands
                @printf(io, "%5d%5d%18.12f\n", ib, ik, eigenvalues[ik][ib])
            end
        end
    end
end

function write_eig(
    filename::AbstractString, eigenvalues::AbstractVector, ::FortranBinaryStream
)
    n_kpts = length(eigenvalues)
    @assert n_kpts > 0 "Empty eigenvalues"
    n_bands = length(eigenvalues[1])

    # gfortran integer is 4 bytes
    Tint = Int32
    # I write in Fortran stream IO, so just plain julia `open`
    open(filename, "w") do io
        for ik in 1:n_kpts
            for ib in 1:n_bands
                write(io, Tint(ib))
                write(io, Tint(ik))
                write(io, Float64(eigenvalues[ik][ib]))
            end
        end
    end
end

function write_eig(filename::AbstractString, eigenvalues::AbstractVector; binary=false)
    if binary
        format = FortranBinaryStream()
    else
        format = FortranText()
    end
    write_eig(filename, eigenvalues, format)

    n_kpts = length(eigenvalues)
    @assert n_kpts > 0 "Empty eigenvalues"
    n_bands = length(eigenvalues[1])
    @info "Writing eig file" filename n_kpts n_bands
end
