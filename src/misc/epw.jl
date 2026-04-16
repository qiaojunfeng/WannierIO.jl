"""
    $(SIGNATURES)

Read the EPW mmn file.

The EPW mmn format is different from that of W90. It does not contain the number
of kpoints/bvectors/bands, so they need to be provided as keyword arguments.

# Arguments
- `file`: The name of the input file, or an `IO`.

# Keyword arguments
- `n_kpts`: number of kpoints
- `n_bvecs`: number of bvectors
- `n_bands`: number of bands

# Return
- `M`: dense array of size `n_bands × n_bands × n_bvecs × n_kpts`
"""
function read_epw_mmn(io::IO; n_kpts::Integer, n_bvecs::Integer, n_bands::Integer)
    M = zeros(ComplexF64, n_bands, n_bands, n_bvecs, n_kpts)

    for ik in 1:n_kpts
        for ib in 1:n_bvecs
            for n in 1:n_bands
                for m in 1:n_bands
                    line = readline(io)
                    line = replace(strip(line), "(" => "", ")" => "", "," => " ")
                    line = split(line)
                    M[m, n, ib, ik] = parse(Float64, line[1]) + parse(Float64, line[2]) * im
                end
            end
        end
    end
    eof(io) ||
        error("Did not reach the end of the file, maybe wrong n_kpts, n_bvecs, or n_bands?")
    return M
end

function read_epw_mmn(
        filename::AbstractString; n_kpts::Integer, n_bvecs::Integer, n_bands::Integer
    )
    return open(filename) do io
        read_epw_mmn(io; n_kpts, n_bvecs, n_bands)
    end
end

"""
Struct for the EPW `.ukk` file.

Similar to the W90 `.chk` file.

$(TYPEDEF)

# Fields

$(FIELDS)
"""
struct Ukk{T <: Real}
    """index of the first band"""
    ibndstart::Int

    """index of the last band"""
    ibndend::Int

    """gauge matrices, size `n_bands * n_wann * n_kpts` array"""
    U::Array{Complex{T}, 3}

    """flag for frozen bands, size `n_bands * n_kpts` matrix"""
    frozen_bands::BitMatrix

    """flag for excluded bands, length-`n_bands + n_excl_bands` vector, where
    `n_excl_bands` is the number of excluded bands"""
    excluded_bands::BitVector

    """centers of WFs, length-`n_wann` vector of `Vec3`.
    Note that EPW uses Cartesian coordinates w.r.t the QE `alat`, so it is dimensionless."""
    centers::Vector{Vec3{T}}
end

n_bands(ukk::Ukk) = size(ukk.U, 1)
n_wannier(ukk::Ukk) = size(ukk.U, 2)
n_kpoints(ukk::Ukk) = size(ukk.U, 3)

"""
    $(SIGNATURES)

Compare two `Ukk` structs.
"""
function Base.isapprox(a::Ukk, b::Ukk)
    return _isapprox(a, b)
end

"""
    $(SIGNATURES)

Read the EPW `.ukk` file.

# Arguments
- `file`: The name of the input file, or an `IO`.

# Return
- `ukk`: the [`Ukk`](@ref) struct
"""
function read_epw_ukk(io::IO)
    lines = readlines(io)

    # Need to 1st read the last part to get the number of WFs
    centers = Vec3{Float64}[]
    for line in Iterators.reverse(lines)
        r = split(strip(line))
        if length(r) == 3
            push!(centers, Vec3(parse.(Float64, r)))
        else
            break
        end
    end
    centers = reverse(centers)
    n_wann = length(centers)
    n_wann > 0 || error("n_wann = $n_wann ≤ 0")

    io2 = IOBuffer(join(lines, "\n") * "\n")
    ibndstart, ibndend = parse.(Int, split(readline(io2)))

    # the unitary matrices
    # now we don't know n_kpts and n_bands yet, so we can only read the
    # complex numbers into a flat vector
    Uflat = ComplexF64[]
    # both the frozen_bands and excluded_bands, still flat vector
    flags = Bool[]

    for line in eachline(io2)
        line = replace(strip(line), "(" => "", ")" => "", "," => " ")
        line = split(line)
        if length(line) == 2
            push!(Uflat, parse(Float64, line[1]) + parse(Float64, line[2]) * im)
        elseif length(line) == 1
            push!(flags, parse_bool(line[1]))
            break
        else
            error("Wrong number of elements in line: $line")
        end
    end

    for line in eachline(io2)
        line = strip(line)
        if length(split(line)) == 1
            push!(flags, parse_bool(line))
        else
            break
        end
    end

    n_kpts_bands_wann = length(Uflat)
    n_kpts_bands = n_kpts_bands_wann ÷ n_wann
    # the kept bands are false in the last part of flags
    excluded_bands = BitVector(flags[(n_kpts_bands + 1):end])
    n_bands = count(!, excluded_bands)
    n_kpts = n_kpts_bands ÷ n_bands
    n_kpts > 0 || error("n_kpts = $n_kpts ≤ 0")
    n_bands > 0 || error("n_bands = $n_bands ≤ 0")

    U = zeros(ComplexF64, n_bands, n_wann, n_kpts)
    counter = 1
    for ik in 1:n_kpts
        for ib in 1:n_bands
            for iw in 1:n_wann
                U[ib, iw, ik] = Uflat[counter]
                counter += 1
            end
        end
    end

    frozen_bands = trues(n_bands, n_kpts)
    counter = 1
    for ik in 1:n_kpts
        for ib in 1:n_bands
            frozen_bands[ib, ik] = flags[counter]
            counter += 1
        end
    end

    return Ukk(
        ibndstart,
        ibndend,
        U,
        frozen_bands,
        excluded_bands,
        centers,
    )
end

function read_epw_ukk(filename::AbstractString)
    return open(filename) do io
        read_epw_ukk(io)
    end
end

"""
    $(SIGNATURES)

Write the EPW `.ukk` file.

# Arguments
- `file`: The name of the output file, or an `IO`.
- `ukk`: the [`Ukk`](@ref) struct

# Examples

See [`Ukk(chk::Chk, alat::Real)`](@ref) for how to construct a `Ukk` from a [`Chk`](@ref).
"""
function write_epw_ukk(io::IO, ukk::Ukk)
    @printf(io, "%d %d\n", ukk.ibndstart, ukk.ibndend)

    # the unitary matrices
    for ik in 1:n_kpoints(ukk)
        for ib in 1:n_bands(ukk)
            for iw in 1:n_wannier(ukk)
                u = ukk.U[ib, iw, ik]
                @printf(io, "(%25.18E,%25.18E)\n", real(u), imag(u))
            end
        end
    end

    # needs also lwindow when disentanglement is used
    for ik in 1:n_kpoints(ukk)
        for ib in 1:n_bands(ukk)
            if ukk.frozen_bands[ib, ik]
                @printf(io, "T\n")
            else
                @printf(io, "F\n")
            end
        end
    end

    for ex in ukk.excluded_bands
        if ex
            @printf(io, "T\n")
        else
            @printf(io, "F\n")
        end
    end

    # now write the Wannier centers to files
    for iw in 1:n_wannier(ukk)
        # meed more precision other WS are not determined properly.
        @printf(io, "%22.12E  %22.12E  %22.12E\n", ukk.centers[iw]...)
    end
    return nothing
end

function write_epw_ukk(filename::AbstractString, ukk::Ukk)
    return open(filename, "w") do io
        write_epw_ukk(io, ukk)
    end
end

"""
    $(SIGNATURES)

Construct a EPW [`Ukk`](@ref) from a W90 [`Chk`](@ref).

# Arguments
- `chk`: the [`Chk`](@ref) struct
- `alat`: the QE `alat` in Å unit. Note that the `alat` from QE stdout file is
    in Bohr unit, you need to do the conversion by multiplying it with
    [`Bohr_QE`](@ref).

# Examples

Convert a W90 `.chk` file to a EPW `.ukk` file:
```julia
using WannierIO
chk = read_chk("BN.chk")
# Note we need QE `alat` for ukk. You can get it
# - either by inspecting the QE stdout file, from line like
#       lattice parameter (alat)  =       6.8330  a.u.
#   where the 6.8330 is the alat in Bohr unit. However, the Bohr constant
#   in W90 and QE are slightly different, to be exact we need to do the unit
#   conversion using QE constant:
alat = 6.8330 * WannierIO.Bohr_QE
# - or better by parsing the QE xml file, and the unit conversion is done automatically
using QuantumEspressoIO: read_pw_xml
alat = read_pw_xml("BN.xml").alat
ukk = Ukk(chk, alat)
WannierIO.write_epw_ukk("BN.ukk", ukk)
```
"""
function Ukk(chk::Chk, alat::Real)
    n_bands = chk.n_bands
    exclude_band_indices = chk.exclude_bands
    n_kpts = chk.n_kpts
    frozen_bands = trues(n_bands, n_kpts)

    n_excl_bands = length(exclude_band_indices)
    n_bands_tot = n_bands + n_excl_bands
    included = trues(n_bands_tot)
    included[exclude_band_indices] .= false
    excluded_bands = .!included
    if n_excl_bands > 0
        ibndstart = findfirst(included)
        ibndend = n_bands_tot - findfirst(reverse(included)) + 1
    else
        ibndstart = 1
        ibndend = n_bands_tot
    end

    # the centers in ukk file is dimensionless: Cartesian coordinates w.r.t alat
    # the centers in chk file is Cartesian coordinates in Å
    # the input arg `alat` should be in Å unit
    centers = chk.r / alat
    Uchk = gauge_matrices(chk)

    return Ukk(
        ibndstart,
        ibndend,
        Uchk,
        frozen_bands,
        excluded_bands,
        centers,
    )
end
