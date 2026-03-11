export read_nnkp, write_nnkp

abstract type Orbital end

"""
Hydrogen-like analytic orbitals.

Follows the definition in Wannier90, see
<https://wannier90.readthedocs.io/en/latest/user_guide/wannier90/postproc/#projections-block>
<https://wannier90.readthedocs.io/en/latest/user_guide/wannier90/projections/>

$(TYPEDEF)

# Fields

$(FIELDS)
"""
@kwdef struct HydrogenOrbital <: Orbital
    """3 real numbers of the projection center, in fractional coordinates"""
    center::Vec3{Float64}
    """positive integer, principle quantum number ``n > 0`` for the radial function"""
    n::Int
    """non-negative integer, angular momentum ``l \\ge 0`` of real spherical
    harmonics ``Y_{lm}(\\theta, \\phi)``"""
    l::Int
    """integer, magnetic quantum number ``m``, ``-l \\leq m \\leq l``"""
    m::Int
    """positive real number, controlling the spread of the radial function"""
    α::Float64
    """3 real numbers, the z-axis from which the polar angle ``\\theta``
    is measured, default is `[0, 0, 1]`"""
    zaxis::Vec3{Float64}
    """3 real numbers, the x-axis from which the azimuthal angle ``\\phi``
    is measured, must be orthogonal to `zaxis`, default is `[1, 0, 0]`"""
    xaxis::Vec3{Float64}
end

Base.NamedTuple(o::HydrogenOrbital) = (; o.center, o.n, o.l, o.m, o.α, o.zaxis, o.xaxis)

"""
    read_nnkp(file)
    read_nnkp(file, ::Wannier90Text)
    read_nnkp(file, ::Wannier90Toml)

Read wannier90 `nnkp` file.

# Arguments
- `file`: The name of the input file, or an `IO`.

# Return
- `lattice`: each column is a lattice vector
- `recip_lattice`: each column is a reciprocal lattice vector
- `kpoints`: length-`n_kpts` vector, each element is `Vec3`, in fractional coordinates
- `projections`: length-`n_projs` vector of `HydrogenOrbital`
- `auto_projections`: optional, the number of Wannier functions `n_wann` for automatic
    initial projections
- `kpb_k`: length-`n_kpts` vector, each element is a length-`n_bvecs` vector of
    integers, index of kpoints
- `kpb_G`: length-`n_kpts` vector, each element is a length-`n_bvecs` vector,
    then each element is `Vec3` for translations, fractional w.r.t `recip_lattice`

Wannier90 `nnkp` file is a plain text format, the 2nd version reads `nnkp` file
in Wannier90 format. The thrid version read a TOML-format `nnkp` file, which is
defined by this package, see [`write_nnkp`](@ref).
The 1st version auto detects the format and parse it.
"""
function read_nnkp end

function read_nnkp(io::IO, ::Wannier90Text)
    read_array(type::Type) = map(x -> parse(type, x), split(readline(io)))

    lattice = zeros(Float64, 3, 3)
    recip_lattice = zeros(Float64, 3, 3)

    n_kpts = nothing
    kpoints = nothing
    projections = nothing
    auto_projections = nothing
    kpb_k = nothing
    kpb_G = nothing

    while !eof(io)
        line = readline(io)

        if occursin("begin real_lattice", line)
            for i in 1:3
                lattice[:, i] = read_array(Float64)
            end

            line = strip(readline(io))
            line == "end real_lattice" || error("`end real_lattice` not found")

        elseif occursin("begin recip_lattice", line)
            for i in 1:3
                recip_lattice[:, i] = read_array(Float64)
            end

            line = strip(readline(io))
            line == "end recip_lattice" || error("`end recip_lattice` not found")

        elseif occursin("begin kpoints", line)
            n_kpts = parse(Int, readline(io))
            kpoints = zeros(Vec3{Float64}, n_kpts)

            for i in 1:n_kpts
                kpoints[i] = Vec3(read_array(Float64))
            end

            line = strip(readline(io))
            line == "end kpoints" || error("`end kpoints` not found")

        elseif occursin("begin projections", line)
            n_projs = parse(Int, readline(io))
            projections = Vector{HydrogenOrbital}(undef, n_projs)
            for i in 1:n_projs
                sline = split(readline(io))
                center = Vec3(parse.(Float64, sline[1:3]))
                l, m, n = parse.(Int, sline[4:6])
                sline = split(readline(io))
                zaxis = Vec3(parse.(Float64, sline[1:3]))
                xaxis = Vec3(parse.(Float64, sline[4:6]))
                α = parse(Float64, sline[7])
                projections[i] = HydrogenOrbital(center, n, l, m, α, zaxis, xaxis)
            end

            line = strip(readline(io))
            line == "end projections" || error("`end projections` not found")

        elseif occursin("begin auto_projections", line)
            # number of WFs
            auto_projections = parse(Int, readline(io))
            # magic number, useless
            i = parse(Int, readline(io))
            i == 0 || error("in `auto_projections`, expected 0, got $i")
            line = strip(readline(io))
            line == "end auto_projections" || error("`end auto_projections` not found")

        elseif occursin("begin nnkpts", line)
            !isnothing(n_kpts) || error("no `kpoints` block before `nnkpts` block?")

            n_bvecs = parse(Int, readline(io))
            kpb_k = [zeros(Int, n_bvecs) for _ in 1:n_kpts]
            kpb_G = [zeros(Vec3{Int}, n_bvecs) for _ in 1:n_kpts]

            for ik in 1:n_kpts
                for ib in 1:n_bvecs
                    arr = read_array(Int)
                    ik == arr[1] || error("expected ik = $ik, got $(arr[1])")
                    kpb_k[ik][ib] = arr[2]
                    kpb_G[ik][ib] = Vec3(arr[3:end])
                end
            end

            line = strip(readline(io))
            line == "end nnkpts" || error("`end nnkpts` not found")
        end
    end

    !iszero(lattice) || error("`real_lattice` not found")
    !iszero(recip_lattice) || error("`recip_lattice` not found")
    !isnothing(kpoints) || error("`kpoints` not found")
    !isnothing(kpb_k) || error("`nnkpts` not found")
    !isnothing(kpb_G) || error("`nnkpts` not found")

    lattice = Mat3(lattice)
    recip_lattice = Mat3(recip_lattice)
    # note I keep the order here: projections first, then auto_projections, ...
    res = OrderedDict{String,Any}()
    isnothing(projections) || (res["projections"] = projections)
    isnothing(auto_projections) || (res["auto_projections"] = auto_projections)
    res["lattice"] = lattice
    res["recip_lattice"] = recip_lattice
    res["kpoints"] = kpoints
    res["kpb_k"] = kpb_k
    res["kpb_G"] = kpb_G
    return res
end

function read_nnkp(io::IO, ::Wannier90Toml)
    # I can just reuse the read_win function, without fix win inputs
    nnkp = read_win(io, Wannier90Toml(); standardize=false)

    # Need to set value to Any, otherwise value type can be too narrow and then
    # I cannot assign Mat3 to it.
    nnkp = OrderedDict{String,Any}(string(k) => v for (k, v) in pairs(nnkp))

    # Some following cleanups
    # Convert to Mat3
    if haskey(nnkp, "lattice") || haskey(nnkp, "recip_lattice")
        for k in ("lattice", "recip_lattice")
            if haskey(nnkp, k)
                nnkp[k] = mat3(nnkp[k])
            end
        end
    end

    # Convert to HydrogenOrbital
    if haskey(nnkp, "projections")
        nnkp["projections"] = map(nnkp["projections"]) do proj
            args = NamedTuple((Symbol(k), v) for (k, v) in proj)
            HydrogenOrbital(; args...)
        end
    end

    # Convert to Vec3
    if haskey(nnkp, "kpb_G")
        nnkp["kpb_G"] = [[Vec3{Int}(G) for G in Gk] for Gk in nnkp["kpb_G"]]
    end

    return nnkp
end

function read_nnkp(filename::AbstractString, format::FileFormat)
    return open(filename) do io
        read_nnkp(io, format)
    end
end

function read_nnkp(file::Union{IO,AbstractString})
    format = istoml(file) ? Wannier90Toml() : Wannier90Text()
    nnkp = read_nnkp(file, format)
    return nnkp
end

"""
    $(SIGNATURES)
"""
@inline function _check_nnkp_required_params(kwargs)
    required_keys = ["lattice", "recip_lattice", "kpoints", "kpb_k", "kpb_G"]
    for k in required_keys
        haskey(kwargs, k) || throw(ArgumentError("Required parameter $k not found"))
    end
end

"""
    write_nnkp(file, params; header)
    write_nnkp(file, params, ::Wannier90Text; header)
    write_nnkp(file, params, ::Wannier90Toml; header)

Write a `nnkp` file that can be used by DFT codes, e.g., QE `pw2wannier90`.

# Arguments
- `file`: The name of the output file, or an `IO`.
- `params`: a `Dict` (or `OrderedDict`) of parameters to be written into the `nnkp` file

The `params` should have at least the following keys:
- `lattice`: each column is a lattice vector
- `recip_lattice`: each column is a reciprocal lattice vector
- `kpoints`: length-`n_kpts` vector of `Vec3`, in fractional coordinates
- `kpb_k`: length-`n_kpts` vector, each element is a length-`n_bvecs` vector of
    integers, index of kpoints
- `kpb_G`: length-`n_kpts` vector, each element is a length-`n_bvecs` vector,
    then each element is a `Vec3` for translation vector, fractional w.r.t. `recip_lattice`

The following keys are optional:
- `projections`: optional, length-`n_projs` vector of `HydrogenOrbital`
- `auto_projections`: optional, the number of Wannier functions `n_wann` for automatic
    initial projections. If given, write an `auto_projections` block
- `exclude_bands`: if given, write the specified band indices in the `exclude_bands` block

# Keyword arguments
- `header`: first line of the file

!!! note

    The variables can also be passed as keyword arguments, e.g.,
    ```julia
    write_nnkp(file; lattice, recip_lattice, kpoints, kpb_k, kpb_G)
    ```
"""
function write_nnkp end

function write_nnkp(io::IO, params::AbstractDict, ::Wannier90Text; header=default_header())
    params = OrderedDict{String,Any}(string(k) => v for (k, v) in pairs(params))
    _check_nnkp_required_params(params)

    lattice = params["lattice"]
    recip_lattice = params["recip_lattice"]
    kpoints = params["kpoints"]
    kpb_k = params["kpb_k"]
    kpb_G = params["kpb_G"]
    projections = get(params, "projections", nothing)
    isnothing(projections) ||
        projections isa AbstractVector{<:HydrogenOrbital} ||
        throw(ArgumentError("projections should be a vector of HydrogenOrbital"))

    auto_projections = get(params, "auto_projections", nothing)
    exclude_bands = get(params, "exclude_bands", nothing)

    n_kpts = length(kpoints)
    n_bvecs = length(kpb_k[1])

    size(recip_lattice) == (3, 3) ||
        throw(DimensionMismatch("size(recip_lattice) != (3, 3)"))
    _check_dimensions_kpb(kpb_k, kpb_G)
    n_kpts == length(kpb_k) ||
        throw(DimensionMismatch("kpoints and kpb_k have different length"))
    size(lattice) == (3, 3) || throw(DimensionMismatch("size(lattice) != (3, 3)"))
    size(recip_lattice) == (3, 3) ||
        throw(DimensionMismatch("size(recip_lattice) != (3, 3)"))

    @printf(io, "%s\n", header)
    @printf(io, "\n")
    # mysterious line, seems not used in W90
    @printf(io, "calc_only_A  :  F\n")
    @printf(io, "\n")

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

    if !isnothing(projections)
        @printf(io, "begin projections\n")
        n_projs = length(projections)
        @printf(io, "%d\n", n_projs)
        for p in projections
            @printf(io, "%12.7f %12.7f %12.7f %3d %3d %3d\n", p.center..., p.l, p.m, p.n)
            @printf(io, "%12.7f %12.7f %12.7f    ", p.zaxis...)
            @printf(io, "%12.7f %12.7f %12.7f    ", p.xaxis...)
            @printf(io, "%8.4f\n", p.α)
        end
        @printf(io, "end projections\n")
        @printf(io, "\n")
    end

    if !isnothing(auto_projections)
        @printf(io, "begin auto_projections\n")
        @printf(io, "%d\n", auto_projections)
        @printf(io, "%d\n", 0)  # magic number, useless
        @printf(io, "end auto_projections\n")
        @printf(io, "\n")
    end

    @printf(io, "begin nnkpts\n")
    @printf(io, "%d\n", n_bvecs)
    for ik in 1:n_kpts
        for ib in 1:n_bvecs
            @printf(io, "%6d %6d %6d %6d %6d\n", ik, kpb_k[ik][ib], kpb_G[ik][ib]...)
        end
    end
    @printf(io, "end nnkpts\n")
    @printf(io, "\n")

    @printf(io, "begin exclude_bands\n")
    if isnothing(exclude_bands)
        @printf(io, "%d\n", 0)
    else
        @printf(io, "%d\n", length(exclude_bands))
        for b in exclude_bands
            @printf(io, "%d\n", b)
        end
    end
    @printf(io, "end exclude_bands\n")
    @printf(io, "\n")
    return nothing
end

function write_nnkp(io::IO, params::AbstractDict, ::Wannier90Toml; header=default_header())
    params = OrderedDict{String,Any}(string(k) => v for (k, v) in pairs(params))
    _check_nnkp_required_params(params)
    kpb_k = params["kpb_k"]
    _check_dimensions_kpb(kpb_k, params["kpb_G"])
    println(io, header, "\n")
    # Note that this requires https://github.com/JuliaLang/julia/pull/57584
    # otherwise it will fail at writing toml file when the `projections`
    # block exists. I.e., this works for julia >= v"1.11.4".
    write_toml(io, params)
    return nothing
end

function write_nnkp(
    filename::AbstractString,
    params::AbstractDict,
    format::FileFormat;
    header=default_header(),
)
    open(filename, "w") do io
        write_nnkp(io, params, format; header)
    end
end

function write_nnkp(
    file::Union{IO,AbstractString}, params::AbstractDict; header=default_header()
)
    format = Wannier90Text()
    write_nnkp(file, params, format; header)
end
