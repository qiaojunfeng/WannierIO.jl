export read_w90_wsvec, write_w90_wsvec

"""
Container for `prefix_wsvec.dat` data.

$(TYPEDEF)

# Fields

$(FIELDS)
"""
struct WsvecDat{IT<:Integer}
    """Whether MDRS interpolation is enabled.
    i.e. the `use_ws_distance` in the header.
    """
    mdrs::Bool

    "The ``\\mathbf{R}`` vectors"
    Rvectors::Vector{Vec3{IT}}

    "MDRS ``\\mathbf{T}_{m n \\mathbf{R}}`` vectors, or `nothing` when `mdrs == false`"
    Tvectors::Union{Vector{Matrix{Vector{Vec3{IT}}}},Nothing}

    "Degeneracies of MDRS ``\\mathbf{T}_{m n \\mathbf{R}}`` vectors, or `nothing` when `mdrs == false`"
    Tdegens::Union{Vector{Matrix{IT}},Nothing}

    "Number of Wannier functions"
    n_wann::Int

    "Header line"
    header::String
end

"""
For Wigner-Seitz Rvectors, needs to provide a `n_wann` for number of Wannier functions.
"""
function WsvecDat(Rvectors::Vector{<:Vec3}, n_wann::Integer, header::String)
    return WsvecDat(false, Rvectors, nothing, nothing, n_wann, header)
end

"""
For MDRS Rvectors, the `n_wann` is optional and can be automatically determined from the `Tvectors`.
"""
function WsvecDat(
    Rvectors::Vector{<:Vec3},
    Tvectors::Vector{<:Matrix},
    Tdegens::Vector{<:Matrix},
    header::String,
)
    n_wann = size(Tvectors[1], 1)
    return WsvecDat(true, Rvectors, Tvectors, Tdegens, n_wann, header)
end

"""
    $(SIGNATURES)

Read `prefix_wsvec.dat`.

# Return
- [`WsvecDat`](@ref) struct containing the data in the file
"""
function read_w90_wsvec(io::IO)
    header = strip(readline(io))

    # check `use_ws_distance`
    mdrs = false
    mdrs_str = split(header)[end]
    if occursin("use_ws_distance=", mdrs_str)
        mdrs_str = lowercase(split(header, "use_ws_distance=")[2])
        mdrs = parse_bool(mdrs_str)
    end

    Rmn = Vector{Vector{Int}}()
    # Tvectors_flat[iRmn][iT] is Vec3 for T-vector, where iRmn is the index for the
    # combination of Rvector and (m, n), iT is the index of T-vector at iRmn.
    # Later on we will reorder this flattened vector, i.e., unfold the
    # iRmn index into (iR, m, n).
    Tvectors_flat = Vector{Vector{Vec3{Int}}}()
    Tdegens_flat = Vector{Int}()

    while !eof(io)
        line = strip(readline(io))
        # the last line is empty
        if length(line) == 0
            continue
        end

        Rx, Ry, Rz, m, n = parse.(Int, split(line))
        push!(Rmn, [Rx, Ry, Rz, m, n])

        n_T = parse(Int, strip(readline(io)))
        push!(Tdegens_flat, n_T)

        T = zeros(Vec3{Int}, n_T)
        for iT in 1:n_T
            line = strip(readline(io))
            T[iT] = Vec3(parse.(Int, split(line))...)
        end
        push!(Tvectors_flat, T)
    end

    # get number of WFs
    n_wann = length(unique(i[end] for i in Rmn))
    n_Rvecs = length(Rmn) ÷ n_wann^2

    Rvectors = zeros(Vec3{Int}, n_Rvecs)
    iR = 1
    for rmn in Rmn
        m, n = rmn[4:5]
        if m == 1 && n == 1
            Rvectors[iR] = Vec3(rmn[1:3])
            iR += 1
        end
    end

    if !mdrs
        return WsvecDat(Rvectors, n_wann, String(header))
    end

    # Objective: reorder Tvectors_flat -> Tvectors, Tdegens_flat -> Tdegens
    # such that Tvectors[iR][m, n][iT] is Vec3 for T-vector
    Tvectors = [Matrix{Vector{Vec3{Int}}}(undef, n_wann, n_wann) for _ in 1:n_Rvecs]
    # and Tdegens[iR][m, n] is the degeneracy of T-vector at iR-th Rvector &
    # between m-th and n-th WFs
    Tdegens = [Matrix{Int}(undef, n_wann, n_wann) for _ in 1:n_Rvecs]
    iR = 1
    for (iRmn, rmn) in enumerate(Rmn)
        m, n = rmn[(end - 1):end]
        Tvectors[iR][m, n] = Tvectors_flat[iRmn]
        Tdegens[iR][m, n] = Tdegens_flat[iRmn]
        # next Rvector
        if m == n_wann && n == n_wann
            iR += 1
        end
    end

    return WsvecDat(Rvectors, Tvectors, Tdegens, String(header))
end

function read_w90_wsvec(filename::AbstractString)
    return open(filename) do io
        read_w90_wsvec(io)
    end
end

"""
    $(SIGNATURES)

Write `prefix_wsvec.dat`.
"""
function write_w90_wsvec(io::IO, wsvec::WsvecDat)
    if wsvec.mdrs
        println(io, wsvec.header * "  with use_ws_distance=.true.")
    else
        println(io, wsvec.header * "  with use_ws_distance=.false.")
    end

    n_Rvecs = length(wsvec.Rvectors)
    for iR in 1:n_Rvecs
        R = wsvec.Rvectors[iR]
        for m in 1:wsvec.n_wann
            for n in 1:wsvec.n_wann
                @printf(io, "%5d %5d %5d %5d %5d\n", R..., m, n)

                if wsvec.mdrs
                    @printf(io, "%5d\n", wsvec.Tdegens[iR][m, n])
                    for T in wsvec.Tvectors[iR][m, n]
                        @printf(io, "%5d %5d %5d\n", T...)
                    end
                else
                    @printf(io, "%5d\n", 1)
                    @printf(io, "%5d %5d %5d\n", 0, 0, 0)
                end
            end
        end
    end

    return nothing
end

function write_w90_wsvec(filename::AbstractString, wsvec::WsvecDat)
    open(filename, "w") do io
        write_w90_wsvec(io, wsvec)
    end
end
