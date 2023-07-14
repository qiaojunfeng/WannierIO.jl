
export read_spn, write_spn

"""
    read_spn(filename)
    read_spn(filename, ::FortranText)

Read the wannier90 `spn` file.

# Return
- `Sx`: spin x, a length-`n_kpts` vector, each element is a `n_bands`-by-`n_bands` matrix
- `Sy`: spin y, a length-`n_kpts` vector, each element is a `n_bands`-by-`n_bands` matrix
- `Sz`: spin z, a length-`n_kpts` vector, each element is a `n_bands`-by-`n_bands` matrix
- `header`: 1st line of the file
"""
function read_spn end

function read_spn(filename::AbstractString, ::FortranText)
    spn = open("$filename") do io
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
        @assert eof(io)

        return (; Sx, Sy, Sz, header)
    end
    return spn
end

function read_spn(filename::AbstractString)
    if isbinary(filename)
        format = FortranBinary()
    else
        format = FortranText()
    end
    Sx, Sy, Sz, header = read_spn(filename, format)

    n_kpts = length(Sx)
    @assert n_kpts > 0 "empty spn matrix"
    n_bands = size(Sx[1], 1)
    @info "Reading spn file" filename header n_kpts n_bands

    return Sx, Sy, Sz
end

"""
    write_spn(filename, Sx, Sy, Sz; binary=false, header)
    write_spn(filename, Sx, Sy, Sz, ::FortranText; header)

Write the `spn` file.
"""
function write_spn end

@inline function _check_dimensions_Sx_Sy_Sz(Sx, Sy, Sz)
    n_kpts = length(Sx)
    @assert n_kpts > 0 "empty spn matrix"
    @assert n_kpts == length(Sy) == length(Sz) "Sx, Sy, Sz must have the same length"
    n_bands = size(Sx[1], 1)
    @assert all(size.(Sx) .== Ref((n_bands, n_bands))) "Sx[ik] must be a square matrix"
    @assert all(size.(Sy) .== Ref((n_bands, n_bands))) "Sy[ik] must be a square matrix"
    @assert all(size.(Sz) .== Ref((n_bands, n_bands))) "Sz[ik] must be a square matrix"
end

function write_spn(
    filename::AbstractString,
    Sx::AbstractVector,
    Sy::AbstractVector,
    Sz::AbstractVector,
    ::FortranText;
    header=default_header(),
)
    _check_dimensions_Sx_Sy_Sz(Sx, Sy, Sz)
    n_kpts = length(Sx)
    n_bands = size(Sx[1], 1)

    open(filename, "w") do io
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
    end
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
    n_kpts = length(Sx)
    n_bands = size(Sx[1], 1)
    @info "Writing spn file" filename header n_kpts n_bands

    return write_spn(filename, Sx, Sy, Sz; header)
end
