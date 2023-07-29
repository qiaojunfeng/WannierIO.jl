export read_w90_wsvec, write_w90_wsvec

"""
    $(SIGNATURES)

Read `prefix_wsvec.dat`.

# Return
- `mdrs`: whether use MDRS interpolation, i.e. the `use_ws_distance` in the header
- `Rvectors`: the ``\\mathbf{R}``-vectors
- `Tvectors`: the ``\\mathbf{T}_{m n \\mathbf{R}}``-vectors.
    Returned only `mdrs = true`.
- `Tdegens`: the degeneracies of ``\\mathbf{T}_{m n \\mathbf{R}}``-vectors.
    Returned only `mdrs = true`.
- `n_wann`: number of WFs
- `header`: the first line of the file
"""
function read_w90_wsvec(filename::AbstractString)
    header, mdrs, Rmn, Tvectors_flat, Tdegens_flat = open(filename) do io
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

        return header, mdrs, Rmn, Tvectors_flat, Tdegens_flat
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
        return (; mdrs, Rvectors, n_wann, header)
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

    @info "Reading wsvec.dat file" filename header mdrs n_wann n_Rvecs
    return (; mdrs, Rvectors, Tvectors, Tdegens, n_wann, header)
end

"""
    $(SIGNATURES)

Write `prefix_wsvec.dat`.

# Keyword Arguments
- `n_wann`: for Wigner-Seitz Rvectors, needs to provide a `n_wann` for number of
    Wannier functions; for MDRS Rvectors, the `n_wann` is optional and can be
    automatically determined from the `Tvectors`
- `Tvectors` and `Tdegens`: if provided, write in MDRS format; otherwise, write in
    Wigner-Seitz format
Also see the return values of [`read_w90_wsvec`](@ref).
"""
function write_w90_wsvec(
    filename::AbstractString;
    Rvectors::AbstractVector,
    n_wann::Union{Integer,Nothing}=nothing,
    Tvectors::Union{AbstractVector,Nothing}=nothing,
    Tdegens::Union{AbstractVector,Nothing}=nothing,
    header=default_header(),
)
    @assert (Tvectors !== nothing && Tdegens !== nothing) ||
        (Tvectors === nothing && Tdegens === nothing) "Tvectors and Tdegens must be both nothing or both not nothing"
    mdrs = Tvectors !== nothing
    n_Rvecs = length(Rvectors)
    @assert n_Rvecs > 0 "empty Rvectors"
    if mdrs
        @assert length(Tvectors) == length(Tdegens) == n_Rvecs "different n_Rvecs in Rvectors, Tvectors, and Tdegens"
        if n_wann === nothing
            n_wann = size(Tvectors[1], 1)
        else
            @assert n_wann == size(Tvectors[1], 1) "different n_wann in Tvectors and n_wann"
        end
    else
        @assert n_wann !== nothing "n_wann must be provided for Wigner-Seitz format"
    end
    @info "Writing wsvec.dat file" filename header mdrs n_wann n_Rvecs

    open(filename, "w") do io
        if mdrs
            println(io, header * "  with use_ws_distance=.true.")
        else
            println(io, header * "  with use_ws_distance=.false.")
        end

        for iR in 1:n_Rvecs
            R = Rvectors[iR]
            for m in 1:n_wann
                for n in 1:n_wann
                    @printf(io, "%5d %5d %5d %5d %5d\n", R..., m, n)

                    if mdrs
                        @printf(io, "%5d\n", Tdegens[iR][m, n])
                        for T in Tvectors[iR][m, n]
                            @printf(io, "%5d %5d %5d\n", T...)
                        end
                    else
                        @printf(io, "%5d\n", 1)
                        @printf(io, "%5d %5d %5d\n", 0, 0, 0)
                    end
                end
            end
        end
    end
end
