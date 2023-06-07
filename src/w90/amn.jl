using Printf: @printf, @sprintf
using Dates: now

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

    A = [zeros(ComplexF64, n_bands, n_wann) for _ in 1:n_kpts]

    while !eof(io)
        line = readline(io)
        arr = split(line)

        m, n, k = parse.(Int64, arr[1:3])

        a = parse(Float64, arr[4]) + im * parse(Float64, arr[5])
        A[k][m, n] = a
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

    A = [zeros(ComplexF64, n_bands, n_wann) for _ in 1:n_kpts]

    while !eof(io)
        m = read(io, Tint)
        n = read(io, Tint)
        k = read(io, Tint)
        r = read(io, Float64)
        i = read(io, Float64)
        A[k][m, n] = r + im * i
    end

    close(io)

    return header, A
end

"""
    read_amn(filename::AbstractString)

Read the `amn` file.

# Return
- `A`: length-`n_kpts` vector, each element is a `n_bands * n_wann` matrix.
"""
function read_amn(filename::AbstractString)
    @info "Reading $filename"

    if isbinary(filename)
        header, A = _read_amn_bin(filename)
    else
        header, A = _read_amn_fmt(filename)
    end
    n_bands, n_wann = size(A[1])
    n_kpts = length(A)

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
function _write_amn_fmt(filename::AbstractString, A::AbstractVector, header::AbstractString)
    n_bands, n_wann = size(A[1])
    n_kpts = length(A)

    io = open(filename, "w")

    header = strip(header)
    write(io, header, "\n")

    @printf(io, "%3d %4d %4d\n", n_bands, n_kpts, n_wann)

    for ik in 1:n_kpts
        for iw in 1:n_wann
            for ib in 1:n_bands
                a = A[ik][ib, iw]
                @printf(io, "%5d %4d %4d  %16.12f  %16.12f\n", ib, iw, ik, real(a), imag(a))
            end
        end
    end

    return close(io)
end

"""
Write A as Fortran unformatted file.
"""
function _write_amn_bin(filename::AbstractString, A::AbstractVector, header::AbstractString)
    n_bands, n_wann = size(A[1])
    n_kpts = length(A)

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
                a = A[ik][ib, iw]
                write(io, Float64(real(a)))
                write(io, Float64(imag(a)))
            end
        end
    end

    return close(io)
end

"""
    write_amn(filename, A, header; binary=false)

Write `amn` file.

# Arguments
- `filename`: output filename
- `A`: a length-`n_kpts` vector, each element is a `n_bands * n_wann` matrix
- `header`: 1st line of the file

# Keyword arguments
- `binary`: write as Fortran unformatted file
"""
function write_amn(
    filename::AbstractString, A::AbstractVector, header::AbstractString; binary::Bool=false
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

"""
    write_amn(filename, A)

Write `amn` file.

Default header is `Created by WannierIO.jl CURRENT_DATE`.

# Arguments
- `filename`: output filename
- `A`: a length-`n_kpts` vector, each element is a `n_bands * n_wann` matrix

# Keyword arguments
- `binary`: write as Fortran unformatted file
"""
function write_amn(filename::AbstractString, A::AbstractVector; binary::Bool=false)
    header = @sprintf "Created by WannierIO.jl %s" string(now())
    return write_amn(filename, A, header; binary=binary)
end
