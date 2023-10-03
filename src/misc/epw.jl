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

    """centers of WFs, length-`n_wann` vector of `Vec3`"""
    centers::Vector{Vec3{T}}
end

"""
    $(SIGNATURES)

Write the EPW `.ukk` file.

# Arguments
- `filename`: the output file name
- `ukk`: the [`Ukk`](@ref) struct
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
"""
function Ukk(chk::Chk)
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

    centers = chk.r
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

"""
    $(SIGNATURES)

Write a EPW `.ukk` file from a W90 [`Chk`](@ref).
"""
function write_epw_ukk(filename::AbstractString, chk::Chk)
    ukk = Ukk(chk)
    return write_epw_ukk(filename, ukk)
end
