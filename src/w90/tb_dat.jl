export read_w90_tb_dat, write_w90_tb_dat

"""
Container for `prefix_tb.dat` data.

$(TYPEDEF)

# Fields

$(FIELDS)
"""
struct TbDat{T <: Real, IT <: Integer}
    "Header line"
    header::String

    "Lattice matrix, columns are lattice vectors in Å"
    lattice::Mat3{T}

    "``\\mathbf{R}``-vectors on which operators are defined"
    Rvectors::Vector{Vec3{IT}}

    "Degeneracies of each ``\\mathbf{R}``-vector"
    Rdegens::Vector{IT}

    "Hamiltonian matrices ``\\mathbf{H}(\\mathbf{R})``"
    H::Vector{Matrix{Complex{T}}}

    "x-component of position operator"
    r_x::Vector{Matrix{Complex{T}}}

    "y-component of position operator"
    r_y::Vector{Matrix{Complex{T}}}

    "z-component of position operator"
    r_z::Vector{Matrix{Complex{T}}}
end

function Base.show(io::IO, tbdat::TbDat)
    n_wann = isempty(tbdat.H) ? 0 : size(tbdat.H[1], 1)
    return print(io, "TbDat(n_Rvecs=$(length(tbdat.H)), n_wann=$(n_wann))")
end

function Base.show(io::IO, ::MIME"text/plain", tbdat::TbDat)
    n_Rvecs = length(tbdat.H)
    n_wann = isempty(tbdat.H) ? 0 : size(tbdat.H[1], 1)
    degen_min = length(tbdat.Rdegens) == 0 ? 0 : minimum(tbdat.Rdegens)
    degen_max = length(tbdat.Rdegens) == 0 ? 0 : maximum(tbdat.Rdegens)

    return print(
        io,
        """TbDat(
          header: $(tbdat.header)
          lattice (Å):
            $(tbdat.lattice[:, 1])
            $(tbdat.lattice[:, 2])
            $(tbdat.lattice[:, 3])
          n_Rvecs: $(n_Rvecs)
          n_wann: $(n_wann)
          Rdegens range: [$(degen_min), ..., $(degen_max)]
          H: Vector{Matrix{$(eltype(tbdat.H[1]) <: Complex ? "Complex" : "Real")}($(n_wann)×$(n_wann))}
          r_x, r_y, r_z: Vector{Matrix{Complex}}($(n_wann)×$(n_wann))
        )""",
    )
end

"""
    $(SIGNATURES)

Read `prefix_tb.dat`.

# Return
- [`TbDat`](@ref) struct containing the data in the file
"""
function read_w90_tb_dat(io::IO)
    header = strip(readline(io))

    # Å unit
    a1 = parse_vector(readstrip(io))
    a2 = parse_vector(readstrip(io))
    a3 = parse_vector(readstrip(io))
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
        line == "" || error("line is not empty")
        Rvectors[iR] = vec3(parse_vector(readstrip(io), Int))

        for n in 1:n_wann
            for m in 1:n_wann
                line = split(readstrip(io))
                m == parse(Int, line[1]) || error(line)
                n == parse(Int, line[2]) || error(line)

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
        line == "" || error("line is not empty")
        Rvectors[iR] == vec3(parse_vector(readstrip(io), Int)) ||
            error("different R vector")
        for n in 1:n_wann
            for m in 1:n_wann
                line = split(readstrip(io))
                m == parse(Int, line[1]) || error("inconsistent m index")
                n == parse(Int, line[2]) || error("inconsistent n index")

                f = parse.(Float64, line[3:8])
                r_x[iR][m, n] = f[1] + im * f[2]
                r_y[iR][m, n] = f[3] + im * f[4]
                r_z[iR][m, n] = f[5] + im * f[6]
            end
        end
    end

    return TbDat(String(header), lattice, Rvectors, Rdegens, H, r_x, r_y, r_z)
end

function read_w90_tb_dat(filename::AbstractString)
    return open(filename) do io
        read_w90_tb_dat(io)
    end
end

"""
    $(SIGNATURES)

Write `prefix_tb.dat`.

# Arguments
See the fields of [`TbDat`](@ref).
"""
function write_w90_tb_dat(io::IO, tbdat::TbDat)
    n_Rvecs = length(tbdat.H)
    n_Rvecs > 0 || error("empty H")
    n_wann = size(tbdat.H[1], 1)

    println(io, strip(tbdat.header))

    # Å unit
    @printf(io, "%21.16f    %21.16f    %21.16f\n", tbdat.lattice[:, 1]...)
    @printf(io, "%21.16f    %21.16f    %21.16f\n", tbdat.lattice[:, 2]...)
    @printf(io, "%21.16f    %21.16f    %21.16f\n", tbdat.lattice[:, 3]...)

    @printf(io, "%d\n", n_wann)
    @printf(io, "%d\n", n_Rvecs)

    n_columns = 15
    for (iR, degen) in enumerate(tbdat.Rdegens)
        @printf(io, "%5d", degen)
        if mod(iR, n_columns) == 0
            println(io)
        end
    end
    if mod(n_Rvecs, n_columns) != 0
        println(io)
    end

    # Hamiltonian
    for iR in 1:n_Rvecs
        @printf(io, "\n%5d %5d %5d\n", tbdat.Rvectors[iR]...)

        for n in 1:n_wann
            for m in 1:n_wann
                reH = real(tbdat.H[iR][m, n])
                imH = imag(tbdat.H[iR][m, n])
                @printf(io, "%5d %5d   %15.8e %15.8e\n", m, n, reH, imH)
            end
        end
    end

    # WF position operator
    for iR in 1:n_Rvecs
        @printf(io, "\n%5d %5d %5d\n", tbdat.Rvectors[iR]...)

        for n in 1:n_wann
            for m in 1:n_wann
                x = tbdat.r_x[iR][m, n]
                y = tbdat.r_y[iR][m, n]
                z = tbdat.r_z[iR][m, n]
                @printf(
                    io,
                    "%5d %5d   %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n",
                    m,
                    n,
                    real(x),
                    imag(x),
                    real(y),
                    imag(y),
                    real(z),
                    imag(z)
                )
            end
        end
    end

    return nothing
end

function write_w90_tb_dat(filename::AbstractString, tbdat::TbDat)
    return open(filename, "w") do io
        write_w90_tb_dat(io, tbdat)
    end
end
