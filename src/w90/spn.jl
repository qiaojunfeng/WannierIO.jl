
export read_spn, write_spn

"""
    read_spn(filename)
    read_spn(filename, ::FortranText)
    read_spn(filename, ::FortranBinary)

Read the wannier90 `spn` file.

# Return
- `Sx`: spin x, a length-`n_kpts` vector, each element is a `n_bands`-by-`n_bands` matrix
- `Sy`: spin y, a length-`n_kpts` vector, each element is a `n_bands`-by-`n_bands` matrix
- `Sz`: spin z, a length-`n_kpts` vector, each element is a `n_bands`-by-`n_bands` matrix
- `header`: 1st line of the file
"""
function read_spn end

function read_spn(io::IO, ::FortranText)
    header = readline(io)

    arr = split(readline(io))
    n_bands, n_kpts = parse.(Int64, arr[1:2])

    Sx = [zeros(ComplexF64, n_bands, n_bands) for _ in 1:n_kpts]
    Sy = [zeros(ComplexF64, n_bands, n_bands) for _ in 1:n_kpts]
    Sz = [zeros(ComplexF64, n_bands, n_bands) for _ in 1:n_kpts]

    for ik in 1:n_kpts
        for m in 1:n_bands
            for n in 1:m
                for Si in (Sx, Sy, Sz)
                    line = split(readline(io))
                    x, y = parse.(Float64, line)
                    Si[ik][n, m] = x + im * y
                    # although each diagonal element of `S` should be real,
                    # actually it has a very small imaginary part,
                    # so we skip the conjugation on the diagonal elements.
                    m == n && continue
                    Si[ik][m, n] = x - im * y
                end
            end
        end
    end
    eof(io) || error(
        "Did not reach the end of the file, maybe the file is corrupted or not in the correct format",
    )

    return (; Sx, Sy, Sz, header)
end

function read_spn(filename::AbstractString, ::FortranText)
    return open("$filename") do io
        read_spn(io, FortranText())
    end
end

function read_spn(io::FortranFile, ::FortranBinary)

    # strip and read line
    header_len = 60
    header = trimstring(read(io, FString{header_len}))

    # gfortran default integer is 4 bytes
    Tint = Int32
    n_bands, n_kpts = read(io, (Tint, 2))

    Sx = [zeros(ComplexF64, n_bands, n_bands) for _ in 1:n_kpts]
    Sy = [zeros(ComplexF64, n_bands, n_bands) for _ in 1:n_kpts]
    Sz = [zeros(ComplexF64, n_bands, n_bands) for _ in 1:n_kpts]

    # upper triangle part, at each kpoint
    spn_tmp = zeros(ComplexF64, 3, n_bands * (n_bands + 1) ÷ 2)

    for ik in 1:n_kpts
        read(io, spn_tmp)
        counter = 0
        for m in 1:n_bands
            for n in 1:m
                counter += 1
                for (i, Si) in enumerate((Sx, Sy, Sz))
                    Si[ik][n, m] = spn_tmp[i, counter]
                    # although each diagonal element of `S` should be real,
                    # actually it has a very small imaginary part,
                    # so we skip the conjugation on the diagonal elements.
                    m == n && continue
                    Si[ik][m, n] = conj(Si[ik][n, m])
                end
            end
        end
    end
    eof(io) || error(
        "Did not reach the end of the file, maybe the file is corrupted or not in the correct format",
    )
    close(io)

    return (; Sx, Sy, Sz, header)
end

function read_spn(filename::AbstractString, ::FortranBinary)
    io = FortranFile(filename)
    return read_spn(io, FortranBinary())
end

function read_spn(filename::AbstractString)
    if isbinary(filename)
        format = FortranBinary()
    else
        format = FortranText()
    end
    Sx, Sy, Sz, header = read_spn(filename, format)

    n_kpts = length(Sx)
    n_kpts > 0 || error("empty spn matrix")
    return Sx, Sy, Sz
end

"""
    write_spn(filename, Sx, Sy, Sz; binary=false, header)
    write_spn(filename, Sx, Sy, Sz, ::FortranText; header)
    write_spn(filename, Sx, Sy, Sz, ::FortranBinary; header)

Write the `spn` file.
"""
function write_spn end

"""
    $(SIGNATURES)
"""
@inline function _check_dimensions_Sx_Sy_Sz(Sx, Sy, Sz)
    n_kpts = length(Sx)
    n_kpts > 0 || throw(ArgumentError("empty spn matrix"))
    n_kpts == length(Sy) == length(Sz) ||
        throw(DimensionMismatch("Sx, Sy, Sz must have the same length"))
    n_bands = size(Sx[1], 1)
    all(size.(Sx) .== Ref((n_bands, n_bands))) ||
        throw(DimensionMismatch("Sx[ik] must be a square matrix"))
    all(size.(Sy) .== Ref((n_bands, n_bands))) ||
        throw(DimensionMismatch("Sy[ik] must be a square matrix"))
    all(size.(Sz) .== Ref((n_bands, n_bands))) ||
        throw(DimensionMismatch("Sz[ik] must be a square matrix"))
end

function write_spn(
    io::IO,
    Sx::AbstractVector,
    Sy::AbstractVector,
    Sz::AbstractVector,
    ::FortranText;
    header=default_header(),
)
    _check_dimensions_Sx_Sy_Sz(Sx, Sy, Sz)
    n_kpts = length(Sx)
    n_bands = size(Sx[1], 1)

    header = strip(header)
    write(io, header, "\n")

    @printf(io, "%3d %4d\n", n_bands, n_kpts)

    for ik in 1:n_kpts
        for m in 1:n_bands
            for n in 1:m
                for Si in (Sx, Sy, Sz)
                    s = Si[ik][n, m]
                    @printf(io, "%26.16e  %26.16e\n", real(s), imag(s))
                end
            end
        end
    end
    return nothing
end

function write_spn(
    filename::AbstractString,
    Sx::AbstractVector,
    Sy::AbstractVector,
    Sz::AbstractVector,
    ::FortranText;
    header=default_header(),
)
    open(filename, "w") do io
        write_spn(io, Sx, Sy, Sz, FortranText(); header)
    end
    return nothing
end

function write_spn(
    io::FortranFile,
    Sx::AbstractVector,
    Sy::AbstractVector,
    Sz::AbstractVector,
    ::FortranBinary;
    header=default_header(),
)
    _check_dimensions_Sx_Sy_Sz(Sx, Sy, Sz)
    n_kpts = length(Sx)
    n_bands = size(Sx[1], 1)

    header_len = 60
    write(io, FString(header_len, string(strip(header))))

    # gfortran default integer is 4 bytes
    Tint = Int32
    write(io, Tint(n_bands), Tint(n_kpts))

    # upper triangle part, at each kpoint
    spn_tmp = zeros(ComplexF64, 3, n_bands * (n_bands + 1) ÷ 2)

    for ik in 1:n_kpts
        counter = 0
        for m in 1:n_bands
            for n in 1:m
                counter += 1
                for (i, Si) in enumerate((Sx, Sy, Sz))
                    spn_tmp[i, counter] = Si[ik][n, m]
                end
            end
        end
        write(io, spn_tmp)
    end
    close(io)
    return nothing
end

function write_spn(
    filename::AbstractString,
    Sx::AbstractVector,
    Sy::AbstractVector,
    Sz::AbstractVector,
    ::FortranBinary;
    header=default_header(),
)
    io = FortranFile(filename, "w")
    return write_spn(io, Sx, Sy, Sz, FortranBinary(); header)
end

function write_spn(
    filename::AbstractString,
    Sx::AbstractVector,
    Sy::AbstractVector,
    Sz::AbstractVector;
    binary=false,
    header=default_header(),
)
    if binary
        format = FortranBinary()
    else
        format = FortranText()
    end
    _check_dimensions_Sx_Sy_Sz(Sx, Sy, Sz)
    return write_spn(filename, Sx, Sy, Sz, format; header)
end
