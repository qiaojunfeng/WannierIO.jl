"""
    $(SIGNATURES)

Read the EPW mmn file.

The EPW mmn format is different from that of W90. It does not contain the number
of kpoints/bvectors/bands, so they need to be provided as keyword arguments.

# Arguments
- `filename`: the mmn file name

# Keyword arguments
- `n_kpts`: number of kpoints
- `n_bvecs`: number of bvectors
- `n_bands`: number of bands

# Return
- `M`: length-`n_kpts` vector of length-`n_bvecs` vector of `n_bands * n_bands` matrices
"""
function read_epw_mmn(
    filename::AbstractString; n_kpts::Integer, n_bvecs::Integer, n_bands::Integer
)
    return open(filename) do io
        M = [[zeros(ComplexF64, n_bands, n_bands) for _ in 1:n_bvecs] for _ in 1:n_kpts]

        for ik in 1:n_kpts
            for ib in 1:n_bvecs
                for n in 1:n_bands
                    for m in 1:n_bands
                        line = readline(io)
                        line = replace(strip(line), "(" => "", ")" => "", "," => " ")
                        line = split(line)
                        M[ik][ib][m, n] =
                            parse(Float64, line[1]) + parse(Float64, line[2]) * im
                    end
                end
            end
        end
        @assert eof(io) "Did not reach the end of the file, maybe wrong n_kpts, n_bvecs, or n_bands?"
        return M
    end
end

"""
Struct for the EPW `.ukk` file.

Similar to the W90 `.chk` file.

$(TYPEDEF)

# Fields

$(FIELDS)
"""
struct Ukk{T<:Real}
    """index of the first band"""
    ibndstart::Int

    """index of the last band"""
    ibndend::Int

    """number of kpoints"""
    n_kpts::Int

    """number of bands"""
    n_bands::Int

    """number of wannier functions"""
    n_wann::Int

    """gauge matrices, length-`n_kpts` vector, each element is a `n_bands * n_wann` matrix"""
    U::Vector{Matrix{Complex{T}}}

    """flag for frozen bands, length-`n_kpts` vector, each element is a length-`n_bands` vector"""
    frozen_bands::Vector{BitVector}

    """flag for excluded bands, length-`n_bands + n_excl_bands` vector, where
    `n_excl_bands` is the number of excluded bands"""
    excluded_bands::BitVector

    """centers of WFs, length-`n_wann` vector of `Vec3`.
    Note that EPW uses Cartesian coordinates w.r.t the QE `alat`, so it is dimensionless."""
    centers::Vector{Vec3{T}}
end

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
- `filename`: the output file name

# Return
- `ukk`: the [`Ukk`](@ref) struct
"""
function read_epw_ukk(filename::AbstractString)
    # Need to 1st read the last part to get the number of WFs
    centers = open(filename) do io
        centers = Vec3{Float64}[]
        # note Julia 1.8 is required for reverse(eachline(io))
        for line in Iterators.reverse(eachline(io))
            r = split(strip(line))
            if length(r) == 3
                push!(centers, Vec3(parse.(Float64, r)))
            else
                break
            end
        end
        return reverse(centers)
    end
    n_wann = length(centers)
    @assert n_wann > 0 "n_wann = $n_wann ≤ 0"

    ibndstart, ibndend, Uflat, flags = open(filename) do io
        ibndstart, ibndend = parse.(Int, split(readline(io)))

        # the unitary matrices
        # now we don't know n_kpts and n_bands yet, so we can only read the
        # complex numbers into a flat vector
        Uflat = ComplexF64[]
        # both the frozen_bands and excluded_bands, still flat vector
        flags = Bool[]

        for line in eachline(io)
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

        for line in eachline(io)
            line = strip(line)
            if length(split(line)) == 1
                push!(flags, parse_bool(line))
            else
                break
            end
        end

        return ibndstart, ibndend, Uflat, flags
    end

    n_kpts_bands_wann = length(Uflat)
    n_kpts_bands = n_kpts_bands_wann ÷ n_wann
    # the kept bands are false in the last part of flags
    excluded_bands = BitVector(flags[(n_kpts_bands + 1):end])
    n_bands = count(!, excluded_bands)
    n_kpts = n_kpts_bands ÷ n_bands
    @assert n_kpts > 0 "n_kpts = $n_kpts ≤ 0"
    @assert n_bands > 0 "n_bands = $n_bands ≤ 0"
    @info "Reading ukk file" filename n_kpts n_bands n_wann

    U = [zeros(ComplexF64, n_bands, n_wann) for _ in 1:n_kpts]
    counter = 1
    for ik in 1:n_kpts
        for ib in 1:n_bands
            for iw in 1:n_wann
                U[ik][ib, iw] = Uflat[counter]
                counter += 1
            end
        end
    end

    frozen_bands = [trues(n_bands) for _ in 1:n_kpts]
    counter = 1
    for ik in 1:n_kpts
        for ib in 1:n_bands
            frozen_bands[ik][ib] = flags[counter]
        end
    end

    return Ukk(
        ibndstart,
        ibndend,
        n_kpts,
        n_bands,
        n_wann,
        U,
        frozen_bands,
        excluded_bands,
        centers,
    )
end

"""
    $(SIGNATURES)

Write the EPW `.ukk` file.

# Arguments
- `filename`: the output file name
- `ukk`: the [`Ukk`](@ref) struct

# Examples

See [`Ukk(chk::Chk, alat::Real)`](@ref) for how to construct a `Ukk` from a [`Chk`](@ref).
"""
function write_epw_ukk(filename::AbstractString, ukk::Ukk)
    open(filename, "w") do io
        @printf(io, "%d %d\n", ukk.ibndstart, ukk.ibndend)

        # the unitary matrices
        for ik in 1:(ukk.n_kpts)
            for ib in 1:(ukk.n_bands)
                for iw in 1:(ukk.n_wann)
                    u = ukk.U[ik][ib, iw]
                    @printf(io, "(%25.18E,%25.18E)\n", real(u), imag(u))
                end
            end
        end

        # needs also lwindow when disentanglement is used
        for ik in 1:(ukk.n_kpts)
            for ib in 1:(ukk.n_bands)
                if ukk.frozen_bands[ik][ib]
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
        for iw in 1:(ukk.n_wann)
            # meed more precision other WS are not determined properly.
            @printf(io, "%22.12E  %22.12E  %22.12E\n", ukk.centers[iw]...)
        end
    end
    @printf("Written to %s\n", filename)
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
alat = read_qe_xml("BN.xml").alat
ukk = Ukk(chk, alat)
WannierIO.write_epw_ukk("BN.ukk", ukk)
```
"""
function Ukk(chk::Chk, alat::Real)
    n_bands = chk.n_bands
    exclude_band_indices = chk.exclude_bands
    n_kpts = chk.n_kpts
    n_wann = chk.n_wann
    frozen_bands = [trues(n_bands) for _ in 1:n_kpts]

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
    Uchk = get_U(chk)

    return Ukk(
        ibndstart,
        ibndend,
        n_kpts,
        n_bands,
        n_wann,
        Uchk,
        frozen_bands,
        excluded_bands,
        centers,
    )
end
