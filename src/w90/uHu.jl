export read_uHu, write_uHu

"""
    read_uHu(filename)
    read_uHu(filename, ::FortranText; transpose_band_indices=true)
    read_uHu(filename, ::FortranBinary; transpose_band_indices=true)

Read the wannier90 `uHu` file.

# Keyword Arguments
- `transpose_band_indices`: QE pw2wannier90.x writes the matrix in a strange
    transposed manner; if reading a QE-generated `uHu` file, this flag should
    be true to restore the band indices order, so that the returned matrix
    has the correct order, i.e.,
    `uHu[ik][ib1, ib2][m, n]` is
    ``\\langle u_{m, k + b_1}| H | u_{n, k + b_2} \\rangle``

# Return
- `uHu`: a length-`n_kpts` vector, each element is a `n_bvecs * n_bvecs` matrix,
    then each element is a  `n_bands * n_bands` matrix
- `header`: 1st line of the file
"""
function read_uHu end

function read_uHu(filename::AbstractString, ::FortranText; transpose_band_indices=true)
    return open("$filename") do io
        header = readline(io)

        arr = split(readline(io))
        n_bands, n_kpts, n_bvecs = parse.(Int64, arr[1:3])

        uHu = map(1:n_kpts) do _
            map(CartesianIndices((n_bvecs, n_bvecs))) do _
                zeros(ComplexF64, n_bands, n_bands)
            end
        end

        for ik in 1:n_kpts
            for ib2 in 1:n_bvecs
                for ib1 in 1:n_bvecs
                    for n in 1:n_bands
                        for m in 1:n_bands
                            line = split(readline(io))
                            x, y = parse.(Float64, line)
                            uHu[ik][ib1, ib2][m, n] = x + im * y
                        end
                    end
                    if transpose_band_indices
                        uHu[ik][ib1, ib2] .= transpose(uHu[ik][ib1, ib2])
                    end
                end
            end
        end
        @assert eof(io)
        return (; uHu, header)
    end
end

function read_uHu(filename::AbstractString, ::FortranBinary; transpose_band_indices=true)
    io = FortranFile(filename)

    # strip and read line
    header_len = 60
    header = trimstring(read(io, FString{header_len}))

    # gfortran default integer is 4 bytes
    Tint = Int32
    n_bands, n_kpts, n_bvecs = read(io, (Tint, 3))

    uHu = map(1:n_kpts) do _
        map(CartesianIndices((n_bvecs, n_bvecs))) do _
            zeros(ComplexF64, n_bands, n_bands)
        end
    end

    # upper triangle part, at each kpoint, b1, and b2
    uHu_tmp = zeros(ComplexF64, n_bands^2)

    for ik in 1:n_kpts
        for ib2 in 1:n_bvecs
            for ib1 in 1:n_bvecs
                read(io, uHu_tmp)
                uHu[ik][ib1, ib2] .= reshape(uHu_tmp, (n_bands, n_bands))
                if transpose_band_indices
                    uHu[ik][ib1, ib2] .= transpose(uHu[ik][ib1, ib2])
                end
            end
        end
    end
    @assert eof(io)
    close(io)
    return (; uHu, header)
end

function read_uHu(filename::AbstractString; kwargs...)
    if isbinary(filename)
        format = FortranBinary()
    else
        format = FortranText()
    end
    uHu, header = read_uHu(filename, format; kwargs...)

    n_kpts = length(uHu)
    @assert n_kpts > 0 "empty uHu matrix"
    n_bvecs = size(uHu[1], 1)
    @assert n_bvecs > 0 "empty uHu matrix"
    n_bands = size(uHu[1][1, 1], 1)
    @info "Reading uHu file" filename header n_kpts n_bvecs n_bands
    return uHu
end

"""
    write_uHu(filename, uHu; binary=false, header)
    write_uHu(filename, uHu, ::FortranText; header)
    write_uHu(filename, uHu, ::FortranBinary; header)

Write the `uHu` file.

# Keyword Arguments
- `transpose_band_indices`: see [`read_uHu`](@ref)
"""
function write_uHu end

function write_uHu(
    filename::AbstractString,
    uHu::AbstractVector,
    ::FortranText;
    header=default_header(),
    transpose_band_indices=true,
)
    n_kpts = length(uHu)
    @assert n_kpts > 0 "empty uHu matrix"
    n_bvecs = size(uHu[1], 1)
    @assert n_bvecs > 0 "empty uHu matrix"
    n_bands = size(uHu[1][1, 1], 1)

    open(filename, "w") do io
        header = strip(header)
        write(io, header, "\n")
        @printf(io, "%3d %4d %4d\n", n_bands, n_kpts, n_bvecs)

        for ik in 1:n_kpts
            for ib2 in 1:n_bvecs
                for ib1 in 1:n_bvecs
                    if transpose_band_indices
                        out_uHu = transpose(uHu[ik][ib1, ib2])
                    else
                        out_uHu = uHu[ik][ib1, ib2]
                    end
                    for n in 1:n_bands
                        for m in 1:n_bands
                            u = out_uHu[m, n]
                            @printf(io, "%20.10e  %20.10e\n", real(u), imag(u))
                        end
                    end
                end
            end
        end
    end
end

function write_uHu(
    filename::AbstractString,
    uHu::AbstractVector,
    ::FortranBinary;
    header=default_header(),
    transpose_band_indices=true,
)
    n_kpts = length(uHu)
    @assert n_kpts > 0 "empty uHu matrix"
    n_bvecs = size(uHu[1], 1)
    @assert n_bvecs > 0 "empty uHu matrix"
    n_bands = size(uHu[1][1, 1], 1)

    io = FortranFile(filename, "w")

    header_len = 60
    write(io, FString(header_len, string(strip(header))))

    # gfortran default integer is 4 bytes
    Tint = Int32
    write(io, Tint(n_bands), Tint(n_kpts), Tint(n_bvecs))

    # upper triangle part, at each kpoint
    uHu_tmp = zeros(ComplexF64, n_bands^2)

    for ik in 1:n_kpts
        for ib2 in 1:n_bvecs
            for ib1 in 1:n_bvecs
                if transpose_band_indices
                    out_uHu = transpose(uHu[ik][ib1, ib2])
                else
                    out_uHu = uHu[ik][ib1, ib2]
                end
                uHu_tmp .= reshape(out_uHu, (n_bands^2,))
                write(io, uHu_tmp)
            end
        end
    end
    close(io)
    return nothing
end

function write_uHu(
    filename::AbstractString,
    uHu::AbstractVector;
    binary=false,
    header=default_header(),
    kwargs...,
)
    if binary
        format = FortranBinary()
    else
        format = FortranText()
    end
    n_kpts = length(uHu)
    @assert n_kpts > 0 "empty uHu matrix"
    n_bvecs = size(uHu[1], 1)
    @assert n_bvecs > 0 "empty uHu matrix"
    n_bands = size(uHu[1][1, 1], 1)
    @info "Writing uHu file" filename header n_kpts n_bands n_bvecs

    return write_uHu(filename, uHu, format; header, kwargs...)
end
