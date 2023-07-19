export read_amn, write_amn

"""
    read_amn(filename)
    read_amn(filename, ::FortranText)
    read_amn(filename, ::FortranBinaryStream)

Read wannier90 `amn` file.

# Return
- `A`: length-`n_kpts` vector, each element is a `n_bands * n_wann` matrix.
- `header`: first line of the file

Note there are three versions of this function: the 1st one is a wrapper
function that automatically detect the format (text or binary) of the file,
and does some additional pretty printing to give user a quick hint of
the dimensions of the A matrix; it internally calls the 2nd or the 3rd version
for actual reading.

Wannier90 only has Fortran text format for `amn`, however I wrote a custom version
of QE pw2wannier90.x that can output Fortran binary format (using Fortran stream
IO) to save some disk space. The 1st function auto detect the file format so it is
transparent to the user.
"""
function read_amn end

function read_amn(filename::AbstractString, ::FortranText)
    res = open("$filename") do io
        header = strip(readline(io))

        line = split(readline(io))
        n_bands, n_kpts, n_wann = parse.(Int64, line[1:3])

        A = [zeros(ComplexF64, n_bands, n_wann) for _ in 1:n_kpts]

        while !eof(io)
            line = split(readline(io))
            m, n, k = parse.(Int64, line[1:3])
            a = parse(Float64, line[4]) + im * parse(Float64, line[5])
            A[k][m, n] = a
        end

        return (; A, header)
    end

    return res
end

function read_amn(filename::AbstractString, ::FortranBinaryStream)
    # I use stream io to write amn, so I should use plain julia `open`
    # io = FortranFile("$filename")
    res = open("$filename") do io
        header_len = 60
        header = read(io, FString{header_len})
        # From FString to String
        header = strip(String(header))

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

        return (; A, header)
    end

    return res
end

function read_amn(filename::AbstractString)
    if isbinary(filename)
        format = FortranBinaryStream()
    else
        format = FortranText()
    end
    A, header = read_amn(filename, format)

    n_kpts = length(A)
    @assert n_kpts > 0 "A is empty"
    n_bands, n_wann = size(A[1])
    @info "Reading amn file" filename header n_kpts n_bands n_wann

    # I do not return header here, since
    # - it is already printed by @info
    # - user can directly use `A = read_amn(filename)` to load it, without
    #   the need to unpack the NamedTuple
    return A
end

"""
    write_amn(filename, A; header=default_header(), binary=false)
    write_amn(filename, A, ::FortranText; header=default_header())
    write_amn(filename, A, ::FortranBinaryStream; header=default_header())

Write wannier90 `amn` file.

# Arguments
- `filename`: output filename
- `A`: a length-`n_kpts` vector, each element is a `n_bands * n_wann` matrix

# Keyword arguments
- `header`: 1st line of the file
- `binary`: write as Fortran unformatted file, which is the Wannier90 default.
    Here the `binary` kwargs is provided for convenience.

Same as [`read_amn`](@ref) there are three versions of this function, the 1st
one is a wrapper function, it calls the 2nd or the 3rd version depending on
the `binary` kwargs.
"""
function write_amn end

function write_amn(
    filename::AbstractString, A::AbstractVector, ::FortranText; header=default_header()
)
    n_kpts = length(A)
    @assert n_kpts > 0 "A is empty"
    n_bands, n_wann = size(A[1])

    open(filename, "w") do io
        write(io, header, "\n")

        @printf(io, "%3d %4d %4d\n", n_bands, n_kpts, n_wann)

        for ik in 1:n_kpts
            for iw in 1:n_wann
                for ib in 1:n_bands
                    a = A[ik][ib, iw]
                    @printf(
                        io, "%5d %4d %4d  %16.12f  %16.12f\n", ib, iw, ik, real(a), imag(a)
                    )
                end
            end
        end
    end

    return nothing
end

function write_amn(
    filename::AbstractString,
    A::AbstractVector,
    ::FortranBinaryStream;
    header=default_header(),
)
    n_kpts = length(A)
    @assert n_kpts > 0 "A is empty"
    n_bands, n_wann = size(A[1])

    # I write in Fortran stream io format.
    open(filename, "w") do io
        # I need to convert to String instead of SubString, for FString
        header = String(header)
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
    end

    return nothing
end

function write_amn(
    filename::AbstractString, A::AbstractVector; header=default_header(), binary=false
)
    n_kpts = length(A)
    @assert n_kpts > 0 "A is empty"
    n_bands, n_wann = size(A[1])

    @info "Writing amn file" filename header n_kpts n_bands n_wann

    if binary
        format = FortranBinaryStream()
    else
        format = FortranText()
    end
    return write_amn(filename, A, format; header)
end
