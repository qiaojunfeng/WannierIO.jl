export read_w90_wsvec_dat, write_w90_wsvec_dat

"""
Container for `prefix_wsvec.dat` data.

$(TYPEDEF)

# Fields

$(FIELDS)
"""
struct WsvecDat{IT <: Integer}
    "Header line"
    header::String

    """Whether MDRS interpolation is enabled.
    i.e. the `use_ws_distance` in the header.
    """
    mdrs::Bool

    "The ``\\mathbf{R}`` vectors"
    Rvectors::Vector{Vec3{IT}}

    """MDRS ``\\mathbf{T}_{m n \\mathbf{R}}`` vectors of size `n_wann × n_wann × n_Rvecs`,
    or `nothing` when `mdrs == false`"""
    Tvectors::Union{Array{Vector{Vec3{IT}}, 3}, Nothing}

    """Degeneracies of MDRS ``\\mathbf{T}_{m n \\mathbf{R}}`` vectors of size `n_wann × n_wann × n_Rvecs`, or `nothing` when `mdrs == false`"""
    Tdegens::Union{Array{IT, 3}, Nothing}

    "Number of Wannier functions"
    n_wann::Int
end

n_wannier(wsvec::WsvecDat) = wsvec.n_wann
n_Rvectors(wsvec::WsvecDat) = length(wsvec.Rvectors)

function Base.show(io::IO, wsvec::WsvecDat)
    return print(
        io,
        "WsvecDat(n_Rvecs=$(n_Rvectors(wsvec)), n_wann=$(n_wannier(wsvec)), mdrs=$(wsvec.mdrs))",
    )
end

function Base.show(io::IO, ::MIME"text/plain", wsvec::WsvecDat)
    return print(
        io,
        """WsvecDat(
          header: $(wsvec.header)
          mdrs: $(wsvec.mdrs)
          n_Rvecs: $(n_Rvectors(wsvec))
          n_wann: $(n_wannier(wsvec))
          Tvectors: $(isnothing(wsvec.Tvectors) ? "nothing" : "present")
          Tdegens: $(isnothing(wsvec.Tdegens) ? "nothing" : "present")
        )""",
    )
end

"""
For Wigner-Seitz Rvectors, needs to provide a `n_wann` for number of Wannier functions.
"""
function WsvecDat(header::AbstractString, Rvectors::AbstractVector{<:AbstractVector}, n_wann::Integer)
    IT = eltype(eltype(Rvectors))
    return WsvecDat{IT}(String(header), false, collect(Vec3{IT}.(Rvectors)), nothing, nothing, n_wann)
end

"""
For MDRS Rvectors, the `n_wann` is optional and can be automatically determined from the `Tvectors`.
"""
function WsvecDat(
        header::AbstractString,
        Rvectors::AbstractVector{<:AbstractVector},
        Tvectors::AbstractArray,
        Tdegens::AbstractArray,
    )
    n_wann = size(Tvectors, 1)
    IT = eltype(eltype(Rvectors))
    return WsvecDat{IT}(
        String(header),
        true,
        collect(Vec3{IT}.(Rvectors)),
        Array{Vector{Vec3{IT}}, 3}(Tvectors),
        Array{IT, 3}(Tdegens),
        n_wann,
    )
end

"""
    $(SIGNATURES)

Read `prefix_wsvec.dat`.

# Return
- [`WsvecDat`](@ref) struct containing the data in the file
"""
function read_w90_wsvec_dat(io::IO)
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
        return WsvecDat(String(header), Rvectors, n_wann)
    end

    # Objective: reorder Tvectors_flat -> Tvectors, Tdegens_flat -> Tdegens
    # such that Tvectors[iR][m, n][iT] is Vec3 for T-vector
    Tvectors = Array{Vector{Vec3{Int}}, 3}(undef, n_wann, n_wann, n_Rvecs)
    # and Tdegens[iR][m, n] is the degeneracy of T-vector at iR-th Rvector &
    # between m-th and n-th WFs
    Tdegens = zeros(Int, n_wann, n_wann, n_Rvecs)
    iR = 1
    for (iRmn, rmn) in enumerate(Rmn)
        m, n = rmn[(end - 1):end]
        Tvectors[m, n, iR] = Tvectors_flat[iRmn]
        Tdegens[m, n, iR] = Tdegens_flat[iRmn]
        # next Rvector
        if m == n_wann && n == n_wann
            iR += 1
        end
    end

    return WsvecDat(String(header), Rvectors, Tvectors, Tdegens)
end

function read_w90_wsvec_dat(filename::AbstractString)
    return open(filename) do io
        read_w90_wsvec_dat(io)
    end
end

"""
    $(SIGNATURES)

Write `prefix_wsvec.dat`.
"""
function write_w90_wsvec_dat(io::IO, wsvec::WsvecDat)
    if wsvec.mdrs
        println(io, wsvec.header * "  with use_ws_distance=.true.")
    else
        println(io, wsvec.header * "  with use_ws_distance=.false.")
    end

    n_Rvecs = n_Rvectors(wsvec)
    n_wann = n_wannier(wsvec)
    for iR in 1:n_Rvecs
        R = wsvec.Rvectors[iR]
        for m in 1:n_wann
            for n in 1:n_wann
                @printf(io, "%5d %5d %5d %5d %5d\n", R..., m, n)

                if wsvec.mdrs
                    @printf(io, "%5d\n", wsvec.Tdegens[m, n, iR])
                    for T in wsvec.Tvectors[m, n, iR]
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

function write_w90_wsvec_dat(filename::AbstractString, wsvec::WsvecDat)
    return open(filename, "w") do io
        write_w90_wsvec_dat(io, wsvec)
    end
end
