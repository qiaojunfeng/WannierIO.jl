using Printf: @printf

export read_eig, write_eig

"""
Reshape a vector of eigenvalues into a matrix of eigenvalues.
"""
function _reshape_eig(
    idx_b::AbstractVector{<:Integer},
    idx_k::AbstractVector{<:Integer},
    eig::AbstractVector{<:Real},
)
    # find unique elements
    n_bands = length(Set(idx_b))
    n_kpts = length(Set(idx_k))
    E = reshape(eig, (n_bands, n_kpts))

    return E
end

"""
Read plain text eig file.
"""
function _read_eig_fmt(filename::AbstractString)
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

    E = _reshape_eig(idx_b, idx_k, eig)
    return E
end

"""
Read binary eig file.
"""
function _read_eig_bin(filename::AbstractString)
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

    E = _reshape_eig(idx_b, idx_k, eig)
    return E
end

"""
    read_eig(filename::AbstractString)

Read the `eig` file.

# Return
- `E`: a `n_bands * n_kpts` array
"""
function read_eig(filename::AbstractString)
    @info "Reading $filename"

    if isbinary(filename)
        E = _read_eig_bin(filename)
    else
        E = _read_eig_fmt(filename)
    end

    n_bands, n_kpts = size(E)

    # check that eigenvalues are in order
    # some times there are small noises
    round_digits(x) = round(x; digits=7)
    for ik in 1:n_kpts
        if !issorted(E[:, ik]; by=round_digits)
            @warn "Eigenvalues are not sorted at " ik E[:, ik]
        end
    end

    println("  n_bands = ", n_bands)
    println("  n_kpts  = ", n_kpts)
    println()

    return E
end

"""
Write plain text eig file.
"""
function _write_eig_fmt(filename::AbstractString, E::AbstractMatrix{<:Real})
    n_bands, n_kpts = size(E)

    open(filename, "w") do io
        for ik in 1:n_kpts
            for ib in 1:n_bands
                @printf(io, "%5d%5d%18.12f\n", ib, ik, E[ib, ik])
            end
        end
    end
end

"""
Write binary eig file.
"""
function _write_eig_bin(filename::AbstractString, E::AbstractMatrix{<:Real})
    n_bands, n_kpts = size(E)

    # gfortran integer is 4 bytes
    Tint = Int32
    # I write in Fortran stream IO, so just plain julia `open`
    open(filename, "w") do io
        for ik in 1:n_kpts
            for ib in 1:n_bands
                write(io, Tint(ib))
                write(io, Tint(ik))
                write(io, Float64(E[ib, ik]))
            end
        end
    end
end

"""
    write_eig(filename::AbstractString, E::AbstractArray; binary=false)

Write `eig` file.

# Arguments
- `E`: `n_bands * n_kpts`

# Keyword arguments
- `binary`: if true write in Fortran binary format
"""
function write_eig(filename::AbstractString, E::AbstractMatrix{<:Real}; binary::Bool=false)
    if binary
        _write_eig_bin(filename, E)
    else
        _write_eig_fmt(filename, E)
    end

    @info "Written to file: $(filename)"
    println()

    return nothing
end
