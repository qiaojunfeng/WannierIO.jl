export read_w90_r_dat, write_w90_r_dat

"""
Container for `prefix_r.dat` data.

$(TYPEDEF)

# Fields

$(FIELDS)
"""
struct RDat{T <: Real, IT <: Integer}
    "Header line"
    header::String

    "``\\mathbf{R}``-vectors on which operators are defined"
    Rvectors::Vector{Vec3{IT}}

    "x-component of position operator, size `n_wann × n_wann × n_Rvecs`"
    rx::Array{Complex{T}, 3}

    "y-component of position operator, size `n_wann × n_wann × n_Rvecs`"
    ry::Array{Complex{T}, 3}

    "z-component of position operator, size `n_wann × n_wann × n_Rvecs`"
    rz::Array{Complex{T}, 3}
end

function RDat(
        header::AbstractString,
        Rvectors::Vector{<:Vec3},
        rx::AbstractVector{<:AbstractMatrix{Complex{T}}},
        ry::AbstractVector{<:AbstractMatrix{Complex{T}}},
        rz::AbstractVector{<:AbstractMatrix{Complex{T}}},
    ) where {T <: Real}
    return RDat(
        String(header),
        collect(Vec3{eltype(first(Rvectors))}.(Rvectors)),
        cat(rx...; dims = 3),
        cat(ry...; dims = 3),
        cat(rz...; dims = 3),
    )
end

function Base.show(io::IO, rdat::RDat)
    n_wann, _, n_Rvecs = size(rdat.rx)
    return print(io, "RDat(n_Rvecs=$(n_Rvecs), n_wann=$(n_wann))")
end

function Base.show(io::IO, ::MIME"text/plain", rdat::RDat)
    n_wann, _, n_Rvecs = size(rdat.rx)

    return print(
        io,
        """RDat(
          header: $(rdat.header)
          n_Rvecs: $(n_Rvecs)
          n_wann: $(n_wann)
                    rx, ry, rz: Array{Complex}($(n_wann)×$(n_wann)×$(n_Rvecs))
        )""",
    )
end

"""
    $(SIGNATURES)

Read `prefix_r.dat`.

# Return
- [`RDat`](@ref) struct containing the data in the file
"""
function read_w90_r_dat(io::IO)
    header = strip(readline(io))
    n_wann = parse(Int, strip(readline(io)))
    n_Rvecs = parse(Int, strip(readline(io)))

    Rvectors = zeros(Vec3{Int}, n_Rvecs)
    rx = zeros(ComplexF64, n_wann, n_wann, n_Rvecs)
    ry = zeros(ComplexF64, n_wann, n_wann, n_Rvecs)
    rz = zeros(ComplexF64, n_wann, n_wann, n_Rvecs)
    for iR in 1:n_Rvecs
        for n in 1:n_wann
            for m in 1:n_wann
                line = split(strip(readline(io)))
                Rvectors[iR] = parse.(Int, line[1:3])
                m == parse(Int, line[4]) || error(line)
                n == parse(Int, line[5]) || error(line)
                rx[m, n, iR] = complex(parse(Float64, line[6]), parse(Float64, line[7]))
                ry[m, n, iR] = complex(parse(Float64, line[8]), parse(Float64, line[9]))
                rz[m, n, iR] = complex(parse(Float64, line[10]), parse(Float64, line[11]))
            end
        end
    end

    return RDat(String(header), Rvectors, rx, ry, rz)
end

function read_w90_r_dat(filename::AbstractString)
    return open(filename) do io
        read_w90_r_dat(io)
    end
end

"""
    $(SIGNATURES)

Write `prefix_r.dat`.

# Arguments
See the fields of [`RDat`](@ref).
"""
function write_w90_r_dat(io::IO, rdat::RDat)
    n_Rvecs = length(rdat.Rvectors)
    n_Rvecs > 0 || throw(ArgumentError("empty Rvectors"))
    n_Rvecs == size(rdat.rx, 3) == size(rdat.ry, 3) == size(rdat.rz, 3) ||
        throw(DimensionMismatch("inconsistent length"))
    n_wann, _, _ = size(rdat.rx)

    println(io, strip(rdat.header))
    @printf(io, "%d\n", n_wann)
    @printf(io, "%d\n", n_Rvecs)

    for iR in 1:n_Rvecs
        for n in 1:n_wann
            for m in 1:n_wann
                @printf(
                    io,
                    " %4d %4d %4d %4d %4d %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f\n",
                    rdat.Rvectors[iR]...,
                    m,
                    n,
                    real(rdat.rx[m, n, iR]),
                    imag(rdat.rx[m, n, iR]),
                    real(rdat.ry[m, n, iR]),
                    imag(rdat.ry[m, n, iR]),
                    real(rdat.rz[m, n, iR]),
                    imag(rdat.rz[m, n, iR]),
                )
            end
        end
    end

    return nothing
end

function write_w90_r_dat(filename::AbstractString, rdat::RDat)
    return open(filename, "w") do io
        write_w90_r_dat(io, rdat)
    end
end
