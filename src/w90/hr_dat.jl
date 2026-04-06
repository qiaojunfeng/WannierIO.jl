export read_w90_hr_dat, write_w90_hr_dat

"""
Container for `prefix_hr.dat` data.

$(TYPEDEF)

# Fields

$(FIELDS)
"""
struct HrDat{T <: Real, IT <: Integer}
    "Header line"
    header::String

    "`R` vectors on which operators are defined"
    Rvectors::Vector{Vec3{IT}}

    "Degeneracies of each `R` vector"
    Rdegens::Vector{IT}

    "Hamiltonian matrices in real space"
    H::Vector{Matrix{Complex{T}}}
end

function Base.show(io::IO, hrdat::HrDat)
    n_wann = isempty(hrdat.H) ? 0 : size(hrdat.H[1], 1)
    return print(io, "HrDat(n_Rvecs=$(length(hrdat.H)), n_wann=$(n_wann))")
end

function Base.show(io::IO, ::MIME"text/plain", hrdat::HrDat)
    n_Rvecs = length(hrdat.H)
    n_wann = isempty(hrdat.H) ? 0 : size(hrdat.H[1], 1)
    degen_min = length(hrdat.Rdegens) == 0 ? 0 : minimum(hrdat.Rdegens)
    degen_max = length(hrdat.Rdegens) == 0 ? 0 : maximum(hrdat.Rdegens)

    return print(
        io,
        """HrDat(
          header: $(hrdat.header)
          n_Rvecs: $(n_Rvecs)
          n_wann: $(n_wann)
          Rdegens range: [$(degen_min), ..., $(degen_max)]
          H: Vector{Matrix{Complex}}($(n_wann)×$(n_wann))
        )""",
    )
end

"""
    $(SIGNATURES)

Read `prefix_hr.dat`.

# Return
- `Rvectors`: ``\\mathbf{R}``-vectors on which operators are defined
- `Rdegens`: degeneracies of each ``\\mathbf{R}``-vector
- `H`: Hamiltonian ``\\mathbf{H}(\\mathbf{R})``
- `header`: the first line of the file
"""
function read_w90_hr_dat(io::IO)
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
                m == parse(Int, line[4]) || error(line)
                n == parse(Int, line[5]) || error(line)
                H[iR][m, n] = complex(parse(Float64, line[6]), parse(Float64, line[7]))
            end
        end
    end

    return HrDat(String(header), Rvectors, Rdegens, H)
end

function read_w90_hr_dat(filename::AbstractString)
    return open(filename) do io
        read_w90_hr_dat(io)
    end
end

"""
    $(SIGNATURES)

Write `prefix_hr.dat`.

# Arguments
See the fields of [`HrDat`](@ref).
"""
function write_w90_hr_dat(io::IO, hrdat::HrDat)
    n_Rvecs = length(hrdat.H)
    n_Rvecs > 0 || throw(ArgumentError("empty H"))
    n_wann = size(hrdat.H[1], 1)

    println(io, strip(hrdat.header))
    @printf(io, "%d\n", n_wann)
    @printf(io, "%d\n", n_Rvecs)

    n_columns = 15
    for (iR, degen) in enumerate(hrdat.Rdegens)
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
                reH = real(hrdat.H[iR][m, n])
                imH = imag(hrdat.H[iR][m, n])
                @printf(
                    io,
                    " %4d %4d %4d %4d %4d %11.6f %11.6f\n",
                    hrdat.Rvectors[iR]...,
                    m,
                    n,
                    reH,
                    imH
                )
            end
        end
    end

    return nothing
end

function write_w90_hr_dat(filename::AbstractString, hrdat::HrDat)
    return open(filename, "w") do io
        write_w90_hr_dat(io, hrdat)
    end
end
