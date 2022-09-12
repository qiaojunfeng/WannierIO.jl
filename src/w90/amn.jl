using Printf: @printf, @sprintf
using Dates: now
using FortranFiles: FortranFile, FString

export read_amn, write_amn

"""
    _read_amn_fmt(filename::AbstractString)

Read Fortran formatted `amn` file.
"""
function _read_amn_fmt(filename::AbstractString)
    io = open("$filename")

    header = readline(io)

    arr = split(readline(io))
    n_bands, n_kpts, n_wann = parse.(Int64, arr[1:3])

    A = zeros(ComplexF64, n_bands, n_wann, n_kpts)

    while !eof(io)
        line = readline(io)
        arr = split(line)

        m, n, k = parse.(Int64, arr[1:3])

        a = parse(Float64, arr[4]) + im * parse(Float64, arr[5])
        A[m, n, k] = a
    end

    close(io)

    return header, A
end

"""
    _read_amn_bin(filename::AbstractString)

Read Fortran binary (unformatted) `amn` file.

Using Fortran stream IO.
"""
function _read_amn_bin(filename::AbstractString)
    # I use stream io to write amn, so I should use plain julia `open`
    # io = FortranFile("$filename")
    io = open("$filename")

    header_len = 60
    header = read(io, FString{header_len})

    # gfortran default integer size = 4
    # https://gcc.gnu.org/onlinedocs/gfortran/KIND-Type-Parameters.html
    Tint = Int32
    n_bands = read(io, Tint)
    n_kpts = read(io, Tint)
    n_wann = read(io, Tint)

    A = zeros(ComplexF64, n_bands, n_wann, n_kpts)

    while !eof(io)
        m = read(io, Tint)
        n = read(io, Tint)
        k = read(io, Tint)
        r = read(io, Float64)
        i = read(io, Float64)
        A[m, n, k] = r + im * i
    end

    close(io)

    return header, A
end

"""
    read_amn(filename::AbstractString)

Read the `amn` file.

# Return
- `A`: a `n_bands * n_wann * n_kpts` array.
"""
function read_amn(filename::AbstractString)
    @info "Reading $filename"

    if isbinary(filename)
        header, A = _read_amn_bin(filename)
    else
        header, A = _read_amn_fmt(filename)
    end
    n_bands, n_wann, n_kpts = size(A)

    println("  header  = ", header)
    println("  n_bands = ", n_bands)
    println("  n_wann  = ", n_wann)
    println("  n_kpts  = ", n_kpts)
    println()

    return A
end

"""
Write A in plain text format.
"""
function _write_amn_fmt(
    filename::AbstractString, A::AbstractArray{<:Complex,3}, header::AbstractString
)
    n_bands, n_wann, n_kpts = size(A)

    io = open(filename, "w")

    header = strip(header)
    write(io, header, "\n")

    @printf(io, "%3d %4d %4d\n", n_bands, n_kpts, n_wann)

    for ik in 1:n_kpts
        for iw in 1:n_wann
            for ib in 1:n_bands
                a = A[ib, iw, ik]
                @printf(io, "%5d %4d %4d  %16.12f  %16.12f\n", ib, iw, ik, real(a), imag(a))
            end
        end
    end

    return close(io)
end

"""
Write A as Fortran unformatted file.
"""
function _write_amn_bin(
    filename::AbstractString, A::AbstractArray{<:Complex,3}, header::AbstractString
)
    n_bands, n_wann, n_kpts = size(A)

    # I write in Fortran stream io format.
    io = open(filename, "w")

    # I need to convert to String instead of SubString, for FString
    header = String(strip(header))
    header_len = 60
    header = FString(header_len, header)
    write(io, header)

    # gfortran default integer is 4 bytes
    Tint = Int32

    write(io, Tint(n_bands))
    write(io, Tint(n_kpts))
    write(io, Tint(n_wann))

    for ik in 1:n_kpts
        for iw in 1:n_wann
            for ib in 1:n_bands
                write(io, Tint(ib))
                write(io, Tint(iw))
                write(io, Tint(ik))
                a = A[ib, iw, ik]
                write(io, Float64(real(a)))
                write(io, Float64(imag(a)))
            end
        end
    end

    return close(io)
end

"""
    write_amn(filename::AbstractString, A::Array{ComplexF64,3})
    write_amn(filename::AbstractString, A::Array{ComplexF64,3}, header::AbstractString)

Output `amn` file.

# Arguments
- `header`: optional, default is "Created by WannierIO.jl CURRENT_DATE"

# Keyword arguments
- `binary`: write as Fortran unformatted file, optional, default is `false`
"""
function write_amn(
    filename::AbstractString,
    A::AbstractArray{<:Complex,3},
    header::AbstractString;
    binary::Bool=false,
)
    if binary
        _write_amn_bin(filename, A, header)
    else
        _write_amn_fmt(filename, A, header)
    end

    @info "Written to file: $(filename)"
    println()

    return nothing
end

function write_amn(
    filename::AbstractString, A::AbstractArray{<:Complex,3}; binary::Bool=false
)
    header = @sprintf "Created by WannierIO.jl %s" string(now())
    return write_amn(filename, A, header; binary=binary)
end
