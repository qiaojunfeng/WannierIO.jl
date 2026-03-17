
export read_spn, write_spn

"""
Container for wannier90 `spn` data.

$(TYPEDEF)

# Fields

$(FIELDS)
"""
struct Spn{T<:Real}
    "Header line"
    header::String

    """Spin x matrices.
    A length-`n_kpts` vector, each element is a `n_bands`-by-`n_bands` matrix.
    """
    Sx::Vector{Matrix{Complex{T}}}

    """Spin y matrices.
    A length-`n_kpts` vector, each element is a `n_bands`-by-`n_bands` matrix.
    """
    Sy::Vector{Matrix{Complex{T}}}

    """Spin z matrices.
    A length-`n_kpts` vector, each element is a `n_bands`-by-`n_bands` matrix.
    """
    Sz::Vector{Matrix{Complex{T}}}
end

function Base.show(io::IO, spn::Spn)
    n_kpts = length(spn.Sx)
    n_bands = n_kpts == 0 ? 0 : size(spn.Sx[1], 1)
    print(io, "Spn(n_kpts=$(n_kpts), n_bands=$(n_bands))")
end

function Base.show(io::IO, ::MIME"text/plain", spn::Spn)
    n_kpts = length(spn.Sx)
    n_bands = n_kpts == 0 ? 0 : size(spn.Sx[1], 1)

    print(
        io,
        """Spn(
          header: $(spn.header)
          n_kpts: $(n_kpts)
          n_bands: $(n_bands)
          Sx, Sy, Sz: Vector{Matrix{Complex}}($(n_bands)×$(n_bands))
        )""",
    )
end

"""
    read_spn(filename)
    read_spn(file, ::FortranText)
    read_spn(file, ::FortranBinary)

Read the wannier90 `spn` file.

# Arguments
- `file`: The name of the input file, or an `IO`.

# Return
- [`Spn`](@ref) struct containing the data in the file
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

    return Spn(String(header), Sx, Sy, Sz)
end

function read_spn(filename::AbstractString, ::FortranText)
    return open(filename) do io
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

    return Spn(String(header), Sx, Sy, Sz)
end

function read_spn(filename::AbstractString, ::FortranBinary)
    io = FortranFile(filename)
    return read_spn(io, FortranBinary())
end

function read_spn(filename::AbstractString)
    format = detect_fortran_format(filename)
    return read_spn(filename, format)
end

"""
    write_spn(filename, spn; binary=false)
    write_spn(file, spn, ::FortranText)
    write_spn(file, spn, ::FortranBinary)

Write the `spn` file.

# Arguments
- `file`: The name of the output file, or an `IO`.
- `spn`: a [`Spn`](@ref) struct
"""
function write_spn end

function write_spn(io::IO, spn::Spn, ::FortranText)
    n_kpts = length(spn.Sx)
    n_kpts > 0 || throw(ArgumentError("empty spn matrix"))
    n_bands = size(spn.Sx[1], 1)

    write(io, strip(spn.header), "\n")

    @printf(io, "%3d %4d\n", n_bands, n_kpts)

    for ik in 1:n_kpts
        for m in 1:n_bands
            for n in 1:m
                for Si in (spn.Sx, spn.Sy, spn.Sz)
                    s = Si[ik][n, m]
                    @printf(io, "%26.16e  %26.16e\n", real(s), imag(s))
                end
            end
        end
    end
    return nothing
end

function write_spn(filename::AbstractString, spn::Spn, ::FortranText)
    open(filename, "w") do io
        write_spn(io, spn, FortranText())
    end
    return nothing
end

function write_spn(io::FortranFile, spn::Spn, ::FortranBinary)
    n_kpts = length(spn.Sx)
    n_kpts > 0 || throw(ArgumentError("empty spn matrix"))
    n_bands = size(spn.Sx[1], 1)

    header_len = 60
    write(io, FString(header_len, string(strip(spn.header))))

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
                for (i, Si) in enumerate((spn.Sx, spn.Sy, spn.Sz))
                    spn_tmp[i, counter] = Si[ik][n, m]
                end
            end
        end
        write(io, spn_tmp)
    end
    close(io)
    return nothing
end

function write_spn(filename::AbstractString, spn::Spn, ::FortranBinary)
    io = FortranFile(filename, "w")
    return write_spn(io, spn, FortranBinary())
end

function write_spn(filename::AbstractString, spn::Spn; binary=false)
    format = fortran_format(; binary)
    return write_spn(filename, spn, format)
end
