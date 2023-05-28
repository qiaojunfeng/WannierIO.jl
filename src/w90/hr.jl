export read_w90_hrdat

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
