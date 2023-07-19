export read_w90_tbdat

"""
    $(SIGNATURES)

Read `prefix_tb.dat`.

# Return
- `lattice`: each column is a lattice vector in Å
- `Rvectors`: ``\\mathbf{R}``-vectors on which operators are defined
- `Rdegens`: degeneracies of each ``\\mathbf{R}``-vector
- `H`: Hamiltonian ``\\mathbf{H}(\\mathbf{R})``
- `r_x`: ``x``-component of position operator
- `r_x`: ``y``-component of position operator
- `r_x`: ``z``-component of position operator
- `header`: the first line of the file
"""
function read_w90_tbdat(filename::AbstractString)
    # convenice function
    ssrline(io) = split(strip(readline(io)))

    return open(filename) do io
        header = strip(readline(io))

        # Å unit
        a1 = parse.(Float64, ssrline(io))
        a2 = parse.(Float64, ssrline(io))
        a3 = parse.(Float64, ssrline(io))
        # column-major
        lattice = Mat3{Float64}(hcat(a1, a2, a3))

        n_wann = parse(Int, strip(readline(io)))
        n_Rvecs = parse(Int, strip(readline(io)))

        Rdegens = parse_vector(io, Int, n_Rvecs)
        Rvectors = zeros(Vec3{Int}, n_Rvecs)
        H = [Matrix{ComplexF64}(undef, n_wann, n_wann) for _ in 1:n_Rvecs]

        # Hamiltonian
        for iR in 1:n_Rvecs
            line = strip(readline(io))  # empty line
            @assert line == ""
            Rvectors[iR] = Vec3(parse.(Int, ssrline(io))...)

            for n in 1:n_wann
                for m in 1:n_wann
                    line = ssrline(io)
                    @assert m == parse(Int, line[1]) line
                    @assert n == parse(Int, line[2]) line

                    reH, imH = parse.(Float64, line[3:4])
                    H[iR][m, n] = reH + im * imH
                end
            end
        end

        # WF position operator
        r_x = [Matrix{ComplexF64}(undef, n_wann, n_wann) for _ in 1:n_Rvecs]
        r_y = [Matrix{ComplexF64}(undef, n_wann, n_wann) for _ in 1:n_Rvecs]
        r_z = [Matrix{ComplexF64}(undef, n_wann, n_wann) for _ in 1:n_Rvecs]

        for iR in 1:n_Rvecs
            line = strip(readline(io))  # empty line
            @assert line == ""
            @assert Rvectors[iR] == Vec3(parse.(Int, ssrline(io))...)
            for n in 1:n_wann
                for m in 1:n_wann
                    line = ssrline(io)
                    @assert m == parse(Int, line[1])
                    @assert n == parse(Int, line[2])

                    f = parse.(Float64, line[3:8])
                    r_x[iR][m, n] = f[1] + im * f[2]
                    r_y[iR][m, n] = f[3] + im * f[4]
                    r_z[iR][m, n] = f[5] + im * f[6]
                end
            end
        end

        @info "Reading tb.dat file" filename header n_wann n_Rvecs
        return (; lattice, Rvectors, Rdegens, H, r_x, r_y, r_z, header)
    end
end
