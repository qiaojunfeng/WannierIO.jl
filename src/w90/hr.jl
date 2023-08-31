export read_w90_hrdat, write_w90_hrdat

"""
    $(SIGNATURES)

Read `prefix_hr.dat`.

# Return
- `Rvectors`: ``\\mathbf{R}``-vectors on which operators are defined
- `Rdegens`: degeneracies of each ``\\mathbf{R}``-vector
- `H`: Hamiltonian ``\\mathbf{H}(\\mathbf{R})``
- `header`: the first line of the file
"""
function read_w90_hrdat(filename::AbstractString)
    return open(filename) do io
        header = strip(readline(io))
        n_wann = parse(Int, strip(readline(io)))
        n_Rvecs = parse(Int, strip(readline(io)))

        Rdegens = parse_vector(io, Int, n_Rvecs)
        Rvectors = zeros(Vec3{Int}, n_Rvecs)
        H = [Matrix{ComplexF64}(undef, n_wann, n_wann) for _ in 1:n_Rvecs]
        for iR in 1:n_Rvecs
            for n in 1:n_wann
                for m in 1:n_wann
                    line = split(strip(readline(io)))
                    Rvectors[iR] = parse.(Int, line[1:3])
                    @assert m == parse(Int, line[4]) line
                    @assert n == parse(Int, line[5]) line
                    H[iR][m, n] = complex(parse(Float64, line[6]), parse(Float64, line[7]))
                end
            end
        end

        @info "Reading hr.dat file" filename header n_wann n_Rvecs
        return (; Rvectors, Rdegens, H, header)
    end
end

"""
    $(SIGNATURES)

Write `prefix_hr.dat`.

# Keyword arguments
See the return values of [`read_w90_hrdat`](@ref).
"""
function write_w90_hrdat(
    filename::AbstractString;
    Rvectors::AbstractVector,
    Rdegens::AbstractVector,
    H::AbstractVector,
    header=default_header(),
)
    n_Rvecs = length(H)
    @assert n_Rvecs > 0 "empty H"
    n_wann = size(H[1], 1)
    @info "Writing hr.dat file" filename header n_wann n_Rvecs

    return open(filename, "w") do io
        println(io, strip(header))
        @printf(io, "%d\n", n_wann)
        @printf(io, "%d\n", n_Rvecs)

        n_columns = 15
        for (iR, degen) in enumerate(Rdegens)
            @printf(io, "%5d", degen)
            if mod(iR, n_columns) == 0
                println(io)
            end
        end
        if mod(n_Rvecs, n_columns) != 0
            println(io)
        end

        for iR in 1:n_Rvecs
            for n in 1:n_wann
                for m in 1:n_wann
                    reH = real(H[iR][m, n])
                    imH = imag(H[iR][m, n])
                    @printf(
                        io,
                        " %4d %4d %4d %4d %4d %11.6f %11.6f\n",
                        Rvectors[iR]...,
                        m,
                        n,
                        reH,
                        imH
                    )
                end
            end
        end
    end
end
