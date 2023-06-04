using Printf: @printf, @sprintf
using Dates: now

export read_nnkp, write_nnkp

"""
    read_nnkp(filename::AbstractString)

Read the `nnkp` file.

# Return
- `recip_lattice`: each column is a reciprocal lattice vector
- `kpoints`: length-`n_kpts` vector of `Vec3`, in fractional coordinates
- `kpb_k`: length-`n_kpts` vector of length-`n_bvecs` vector of integers, index of kpoints
- `kpb_b`: length-`n_kpts` vector of length-`n_bvecs` vector of `Vec3` for bvectors, fractional w.r.t `recip_lattice`
"""
function read_nnkp(filename::AbstractString)
    @info "Reading nnkp file: $filename"

    io = open(filename)

    read_array(type::Type) = map(x -> parse(type, x), split(readline(io)))

    real_lattice = zeros(Float64, 3, 3)
    recip_lattice = similar(real_lattice)

    n_kpts = nothing
    n_bvecs = nothing
    kpoints = nothing
    kpb_k = nothing
    kpb_b = nothing

    while !eof(io)
        line = readline(io)

        if occursin("begin real_lattice", line)
            for i in 1:3
                real_lattice[:, i] = read_array(Float64)
            end

            line = strip(readline(io))
            line != "end real_lattice" && error("expected end real_lattice")

        elseif occursin("begin recip_lattice", line)
            for i in 1:3
                recip_lattice[:, i] = read_array(Float64)
            end

            line = strip(readline(io))
            line != "end recip_lattice" && error("expected end recip_lattice")

        elseif occursin("begin kpoints", line)
            n_kpts = parse(Int, readline(io))
            kpoints = zeros(Vec3{Float64}, n_kpts)

            for i in 1:n_kpts
                kpoints[i] = Vec3(read_array(Float64)...)
            end

            line = strip(readline(io))
            line != "end kpoints" && error("expected end kpoints")

        elseif occursin("begin nnkpts", line)
            n_kpts === nothing && error("no kpoints block before nnkpts block?")

            n_bvecs = parse(Int, readline(io))
            kpb_k = [zeros(Int, n_bvecs) for _ in 1:n_kpts]
            kpb_b = [zeros(Vec3{Int}, n_bvecs) for _ in 1:n_kpts]

            for ik in 1:n_kpts
                for ib in 1:n_bvecs
                    arr = read_array(Int)
                    ik != arr[1] && error("expected ik = $ik, got $(arr[1])")
                    kpb_k[ik][ib] = arr[2]
                    kpb_b[ik][ib] = Vec3(arr[3:end]...)
                end
            end

            line = strip(readline(io))
            line != "end nnkpts" && error("expected end nnkpts")
        end
    end
    close(io)

    println("  n_kpts  = ", n_kpts)
    println("  n_bvecs = ", n_bvecs)
    println()

    real_lattice = Mat3(real_lattice)
    recip_lattice = Mat3(recip_lattice)

    return (; recip_lattice, kpoints, kpb_k, kpb_b)
end

"""
    write_nnkp(filename, bvectors::BVectors, n_wann)

Write a `nnkp` file that can be used by DFT codes, e.g., QE `pw2wannier90`.

# Arguments
- `recip_lattice`: each column is a reciprocal lattice vector
- `kpoints`: length-`n_kpts` vector of `Vec3`, in fractional coordinates
- `kpb_k`: length-`n_kpts` vector of length-`n_bvecs` vector of integers, index of kpoints
- `kpb_b`: length-`n_kpts` vector of length-`n_bvecs` vector of `Vec3` for bvectors, fractional w.r.t `recip_lattice`
- `n_wann`: if given, write a `auto_projections` block
- `exclude_bands`: if given, write the specified band indexes in the `exclude_bands` block

!!! note

    Only a preliminary version, use `auto_projections`.
"""
function write_nnkp(
    filename::AbstractString,
    recip_lattice::AbstractMatrix,
    kpoints::AbstractVector,
    kpb_k::AbstractVector,
    kpb_b::AbstractVector,
    n_wann::Union{Nothing,Integer}=nothing,
    exclude_bands::Union{Nothing,AbstractVector}=nothing,
)
    size(recip_lattice) == (3, 3) || error("size(recip_lattice) != (3, 3)")
    n_kpts = length(kpoints)
    n_bvecs = length(kpb_k[1])

    (length(kpb_k), length(kpb_k[1])) == (n_kpts, n_bvecs) ||
        error("size(kpb_k) != (n_bvecs, n_kpts)")
    (length(kpb_b), length(kpb_b[1])) == (n_kpts, n_bvecs) ||
        error("size(kpb_b) != (n_bvecs, n_kpts)")

    @info "Writing nnkp file: $filename"

    io = open(filename, "w")

    @printf(io, "File written on %s\n", string(now()))
    @printf(io, "\n")
    # mysterious line, seems not used in W90
    @printf(io, "calc_only_A  :  F\n")
    @printf(io, "\n")

    lattice = get_lattice(recip_lattice)
    @printf(io, "begin real_lattice\n")
    for v in eachcol(lattice)
        @printf(io, "%12.7f %12.7f %12.7f\n", v...)
    end
    @printf(io, "end real_lattice\n")
    @printf(io, "\n")

    @printf(io, "begin recip_lattice\n")
    for v in eachcol(recip_lattice)
        @printf(io, "%12.7f %12.7f %12.7f\n", v...)
    end
    @printf(io, "end recip_lattice\n")
    @printf(io, "\n")

    @printf(io, "begin kpoints\n")
    @printf(io, "%d\n", n_kpts)
    for kpt in kpoints
        @printf(io, "%14.8f %14.8f %14.8f\n", kpt...)
    end
    @printf(io, "end kpoints\n")
    @printf(io, "\n")

    if !isnothing(n_wann)
        @printf(io, "begin projections\n")
        @printf(io, "%d\n", 0)
        @printf(io, "end projections\n")
        @printf(io, "\n")

        @printf(io, "begin auto_projections\n")
        @printf(io, "%d\n", n_wann)
        @printf(io, "%d\n", 0)
        @printf(io, "end auto_projections\n")
        @printf(io, "\n")
    end

    @printf(io, "begin nnkpts\n")
    @printf(io, "%d\n", n_bvecs)
    for ik in 1:n_kpts
        for ib in 1:n_bvecs
            @printf(io, "%6d %6d %6d %6d %6d\n", ik, kpb_k[ik][ib], kpb_b[ik][ib]...)
        end
    end
    @printf(io, "end nnkpts\n")
    @printf(io, "\n")

    @printf(io, "begin exclude_bands\n")
    if exclude_bands === nothing
        @printf(io, "%d\n", 0)
    else
        @printf(io, "%d\n", length(exclude_bands))
        for b in exclude_bands
            @printf(io, "%d\n", b)
        end
    end
    @printf(io, "end exclude_bands\n")
    @printf(io, "\n")

    return close(io)
end

"""
Read a nnkp file in TOML format.

Used only in tests.
"""
function _read_nnkp_toml(filename::AbstractString)
    nnkp = _read_win_toml(filename; fix_inputs=false)

    # Convert to Mat3
    if haskey(nnkp, :recip_lattice)
        # Need to set value to Any, otherwise it is Dict{Symbol,Vector},
        # then I cannot assign Mat3 to it.
        nnkp = Dict{Symbol,Any}(pairs(nnkp))
        nnkp[:recip_lattice] = Mat3(nnkp[:recip_lattice])
        nnkp = NamedTuple(nnkp)
    end

    return nnkp
end

"""
Write a nnkp file in TOML format.

Used only in tests.
Internally calls [`_write_win_toml`](@ref).
"""
function _write_nnkp_toml(filename::AbstractString; kwargs...)
    return _write_win_toml(filename; kwargs...)
end
