export read_w90_rdat, write_w90_rdat

"""
    $(SIGNATURES)

Read `prefix_r.dat`.

# Return
- `Rvectors`: ``\\mathbf{R}``-vectors on which operators are defined
- `r_x`: ``x``-component of position operator
- `r_y`: ``y``-component of position operator
- `r_z`: ``z``-component of position operator
- `header`: the first line of the file
"""
function read_w90_rdat(filename::AbstractString)
    return open(filename) do io
        header = strip(readline(io))
        n_wann = parse(Int, strip(readline(io)))
        n_Rvecs = parse(Int, strip(readline(io)))

        Rvectors = zeros(Vec3{Int}, n_Rvecs)
        r_x = [Matrix{ComplexF64}(undef, n_wann, n_wann) for _ in 1:n_Rvecs]
        r_y = [Matrix{ComplexF64}(undef, n_wann, n_wann) for _ in 1:n_Rvecs]
        r_z = [Matrix{ComplexF64}(undef, n_wann, n_wann) for _ in 1:n_Rvecs]
        for iR in 1:n_Rvecs
            for n in 1:n_wann
                for m in 1:n_wann
                    line = split(strip(readline(io)))
                    Rvectors[iR] = parse.(Int, line[1:3])
                    @assert m == parse(Int, line[4]) line
                    @assert n == parse(Int, line[5]) line
                    r_x[iR][m, n] = complex(
                        parse(Float64, line[6]), parse(Float64, line[7])
                    )
                    r_y[iR][m, n] = complex(
                        parse(Float64, line[8]), parse(Float64, line[9])
                    )
                    r_z[iR][m, n] = complex(
                        parse(Float64, line[10]), parse(Float64, line[11])
                    )
                end
            end
        end

        @info "Reading r.dat file" filename header n_wann n_Rvecs
        return (; Rvectors, r_x, r_y, r_z, header)
    end
end

"""
    $(SIGNATURES)

Write `prefix_r.dat`.

# Keyword arguments
See the return values of [`read_w90_rdat`](@ref).
"""
function write_w90_rdat(
    filename::AbstractString;
    Rvectors::AbstractVector,
    r_x::AbstractVector,
    r_y::AbstractVector,
    r_z::AbstractVector,
    header=default_header(),
)
    n_Rvecs = length(Rvectors)
    @assert n_Rvecs > 0 "empty Rvectors"
    @assert n_Rvecs == length(r_x) == length(r_y) == length(r_z) "inconsistent length"
    n_wann = size(r_x[1], 1)
    @info "Writing r.dat file" filename header n_wann n_Rvecs

    return open(filename, "w") do io
        println(io, strip(header))
        @printf(io, "%d\n", n_wann)
        @printf(io, "%d\n", n_Rvecs)

        for iR in 1:n_Rvecs
            for n in 1:n_wann
                for m in 1:n_wann
                    @printf(
                        io,
                        " %4d %4d %4d %4d %4d %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f\n",
                        Rvectors[iR]...,
                        m,
                        n,
                        real(r_x[iR][m, n]),
                        imag(r_x[iR][m, n]),
                        real(r_y[iR][m, n]),
                        imag(r_y[iR][m, n]),
                        real(r_z[iR][m, n]),
                        imag(r_z[iR][m, n]),
                    )
                end
            end
        end
    end
end
