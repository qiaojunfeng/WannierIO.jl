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
    read_nnkp(file, ::W90InputText)
    read_nnkp(file, ::W90InputToml)

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

function read_nnkp(io::IO, ::W90InputText)
    params = OrderedDict{String, Any}()

    while !eof(io)
        line = strip(readline(io))
        isempty(line) && continue

        isblock, block_name = _nnkp_check_line(line)
        !isblock && continue
        _nnkp_parse_block!(params, io, block_name)
    end

    _nnkp_check_required_blocks(params)
    return params
end

"""Read one line from an nnkp block and strip surrounding whitespace."""
function _nnkp_block_nextline(io::IO, block_name::AbstractString)
    eof(io) && error("Error parsing block `$block_name`: unexpected end of file")
    return strip(readline(io))
end

"""Check if a line marks the end of a named block in an nnkp file."""
@inline function _nnkp_block_isend(line::AbstractString, block_name::AbstractString)
    return line == "end $block_name"
end

"""Assert that a line marks the end of a named block."""
@inline function _nnkp_block_mustend(line::AbstractString, block_name::AbstractString)
    return _nnkp_block_isend(line, block_name) || error("`end $block_name` not found")
end

"""Parse a numeric line from nnkp blocks."""
@inline function _nnkp_read_array(line::AbstractString, type::Type)
    return parse.(type, split(line))
end

"""Check whether an nnkp line starts a block and return `(isblock, block_name)`."""
function _nnkp_check_line(line::AbstractString)
    if startswith(lowercase(line), "begin ")
        block_name = strip(lowercase(line[7:end]))
        isempty(block_name) && error("Invalid block line: $line")
        # Force to String to be type stable
        return true, string(block_name)
    end
    return false, string(line)
end

"""Parse `real_lattice` block and return a 3x3 lattice matrix."""
function _nnkp_parse_block_real_lattice(io::IO)
    block_name = "real_lattice"
    lattice = zeros(Float64, 3, 3)
    for i in 1:3
        lattice[:, i] = _nnkp_read_array(_nnkp_block_nextline(io, block_name), Float64)
    end
    _nnkp_block_mustend(_nnkp_block_nextline(io, block_name), block_name)
    return mat3(lattice)
end

"""Parse `recip_lattice` block and return a 3x3 reciprocal lattice matrix."""
function _nnkp_parse_block_recip_lattice(io::IO)
    block_name = "recip_lattice"
    recip_lattice = zeros(Float64, 3, 3)
    for i in 1:3
        recip_lattice[:, i] = _nnkp_read_array(_nnkp_block_nextline(io, block_name), Float64)
    end
    _nnkp_block_mustend(_nnkp_block_nextline(io, block_name), block_name)
    return mat3(recip_lattice)
end

"""Parse `kpoints` block and return a vector of k-points."""
function _nnkp_parse_block_kpoints(io::IO)
    block_name = "kpoints"
    n_kpts = parse(Int, _nnkp_block_nextline(io, block_name))
    kpoints = zeros(Vec3{Float64}, n_kpts)
    for i in eachindex(kpoints)
        kpoints[i] = Vec3(_nnkp_read_array(_nnkp_block_nextline(io, block_name), Float64))
    end
    _nnkp_block_mustend(_nnkp_block_nextline(io, block_name), block_name)
    return kpoints
end

"""Parse `projections` block and return a vector of `HydrogenOrbital`."""
function _nnkp_parse_block_projections(io::IO)
    block_name = "projections"
    n_projs = parse(Int, _nnkp_block_nextline(io, block_name))
    projections = Vector{HydrogenOrbital}(undef, n_projs)
    for i in eachindex(projections)
        sline = split(_nnkp_block_nextline(io, block_name))
        center = Vec3(parse.(Float64, sline[1:3]))
        l, m, n = parse.(Int, sline[4:6])

        sline = split(_nnkp_block_nextline(io, block_name))
        zaxis = Vec3(parse.(Float64, sline[1:3]))
        xaxis = Vec3(parse.(Float64, sline[4:6]))
        α = parse(Float64, sline[7])
        projections[i] = HydrogenOrbital(center, n, l, m, α, zaxis, xaxis)
    end
    _nnkp_block_mustend(_nnkp_block_nextline(io, block_name), block_name)
    return projections
end

"""Parse `auto_projections` block and return the number of Wannier functions."""
function _nnkp_parse_block_auto_projections(io::IO)
    block_name = "auto_projections"
    auto_projections = parse(Int, _nnkp_block_nextline(io, block_name))
    # magic number, useless
    i = parse(Int, _nnkp_block_nextline(io, block_name))
    i == 0 || error("in `auto_projections`, expected 0, got $i")
    _nnkp_block_mustend(_nnkp_block_nextline(io, block_name), block_name)
    return auto_projections
end

"""Parse `nnkpts` block and return `(kpb_k, kpb_G)`."""
function _nnkp_parse_block_nnkpts(io::IO, n_kpts::Int)
    block_name = "nnkpts"
    n_bvecs = parse(Int, _nnkp_block_nextline(io, block_name))

    kpb_k = [zeros(Int, n_bvecs) for _ in 1:n_kpts]
    kpb_G = [zeros(Vec3{Int}, n_bvecs) for _ in 1:n_kpts]

    for ik in 1:n_kpts
        for ib in 1:n_bvecs
            arr = _nnkp_read_array(_nnkp_block_nextline(io, block_name), Int)
            ik == arr[1] || error("expected ik = $ik, got $(arr[1])")
            kpb_k[ik][ib] = arr[2]
            kpb_G[ik][ib] = Vec3(arr[3:end])
        end
    end

    _nnkp_block_mustend(_nnkp_block_nextline(io, block_name), block_name)
    return kpb_k, kpb_G
end

"""Skip unknown nnkp blocks by consuming lines up to the matching `end <block>`."""
function _nnkp_skip_block(io::IO, block_name::AbstractString)
    line = _nnkp_block_nextline(io, block_name)
    while !_nnkp_block_isend(line, block_name)
        line = _nnkp_block_nextline(io, block_name)
    end
    return nothing
end

"""Dispatch to nnkp block parsers and store parsed values in `params`."""
function _nnkp_parse_block!(params::AbstractDict, io::IO, block_name::AbstractString)
    if block_name == "real_lattice"
        params["lattice"] = _nnkp_parse_block_real_lattice(io)
    elseif block_name == "recip_lattice"
        params["recip_lattice"] = _nnkp_parse_block_recip_lattice(io)
    elseif block_name == "kpoints"
        params["kpoints"] = _nnkp_parse_block_kpoints(io)
    elseif block_name == "projections"
        params["projections"] = _nnkp_parse_block_projections(io)
    elseif block_name == "auto_projections"
        params["auto_projections"] = _nnkp_parse_block_auto_projections(io)
    elseif block_name == "nnkpts"
        haskey(params, "kpoints") || error("no `kpoints` block before `nnkpts` block?")
        kpb_k, kpb_G = _nnkp_parse_block_nnkpts(io, length(params["kpoints"]))
        params["kpb_k"] = kpb_k
        params["kpb_G"] = kpb_G
    else
        _nnkp_skip_block(io, block_name)
    end
    return nothing
end

"""Validate required nnkp blocks after parsing."""
function _nnkp_check_required_blocks(params::AbstractDict)
    haskey(params, "lattice") || error("`real_lattice` not found")
    haskey(params, "recip_lattice") || error("`recip_lattice` not found")
    haskey(params, "kpoints") || error("`kpoints` not found")
    haskey(params, "kpb_k") || error("`nnkpts` not found")
    haskey(params, "kpb_G") || error("`nnkpts` not found")

    !iszero(params["lattice"]) || error("`real_lattice` not found")
    !iszero(params["recip_lattice"]) || error("`recip_lattice` not found")
    return nothing
end

function read_nnkp(io::IO, ::W90InputToml)
    nnkp = read_toml(io)

    # Some following cleanups
    # Convert to HydrogenOrbital
    if haskey(nnkp, "projections")
        nnkp["projections"] = map(nnkp["projections"]) do proj
            args = NamedTuple((Symbol(k), v) for (k, v) in proj)
            HydrogenOrbital(; args...)
        end
    end

    return nnkp
end

function read_nnkp(filename::AbstractString, format::AbstractFileFormat)
    return open(filename) do io
        read_nnkp(io, format)
    end
end

function read_nnkp(file::Union{IO, AbstractString})
    format = detect_w90input_format(file)
    nnkp = read_nnkp(file, format)
    return nnkp
end

"""
    $(SIGNATURES)
"""
@inline function _nnkp_check_required_params(kwargs)
    required_keys = ["lattice", "recip_lattice", "kpoints", "kpb_k", "kpb_G"]
    for k in required_keys
        haskey(kwargs, k) || throw(ArgumentError("Required parameter $k not found"))
    end
    return
end

"""
    write_nnkp(file, params; header)
    write_nnkp(file, params, ::W90InputText; header)
    write_nnkp(file, params, ::W90InputToml; header)

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
"""
function write_nnkp end

function write_nnkp(io::IO, params::AbstractDict, ::W90InputText; header = default_header())
    _nnkp_validate_write_params(params)

    _nnkp_write_header(io, header)
    _nnkp_write_block_real_lattice(io, params["lattice"])
    _nnkp_write_block_recip_lattice(io, params["recip_lattice"])
    _nnkp_write_block_kpoints(io, params["kpoints"])

    projections = get(params, "projections", nothing)
    isnothing(projections) || _nnkp_write_block_projections(io, projections)

    auto_projections = get(params, "auto_projections", nothing)
    isnothing(auto_projections) || _nnkp_write_block_auto_projections(io, auto_projections)

    _nnkp_write_block_nnkpts(io, params["kpb_k"], params["kpb_G"])
    _nnkp_write_block_exclude_bands(io, get(params, "exclude_bands", nothing))

    return nothing
end

"""Validate nnkp parameters required by the text writer."""
function _nnkp_validate_write_params(params::AbstractDict)
    _nnkp_check_required_params(params)

    projections = get(params, "projections", nothing)
    isnothing(projections) ||
        projections isa AbstractVector{<:HydrogenOrbital} ||
        throw(ArgumentError("projections should be a vector of HydrogenOrbital"))

    lattice = params["lattice"]
    recip_lattice = params["recip_lattice"]
    kpoints = params["kpoints"]
    kpb_k = params["kpb_k"]
    kpb_G = params["kpb_G"]

    _check_dimensions_kpb(kpb_k, kpb_G)
    length(kpoints) == length(kpb_k) ||
        throw(DimensionMismatch("kpoints and kpb_k have different length"))
    size(lattice) == (3, 3) || throw(DimensionMismatch("size(lattice) != (3, 3)"))
    size(recip_lattice) == (3, 3) ||
        throw(DimensionMismatch("size(recip_lattice) != (3, 3)"))
    return nothing
end

"""Write nnkp header and fixed metadata line."""
function _nnkp_write_header(io::IO, header)
    @printf(io, "%s\n", header)
    @printf(io, "\n")
    # mysterious line, seems not used in W90
    @printf(io, "calc_only_A  :  F\n")
    @printf(io, "\n")
    return nothing
end

"""Write `real_lattice` block."""
function _nnkp_write_block_real_lattice(io::IO, lattice)
    @printf(io, "begin real_lattice\n")
    for v in eachcol(lattice)
        @printf(io, "%12.7f %12.7f %12.7f\n", v...)
    end
    @printf(io, "end real_lattice\n")
    @printf(io, "\n")
    return nothing
end

"""Write `recip_lattice` block."""
function _nnkp_write_block_recip_lattice(io::IO, recip_lattice)
    @printf(io, "begin recip_lattice\n")
    for v in eachcol(recip_lattice)
        @printf(io, "%12.7f %12.7f %12.7f\n", v...)
    end
    @printf(io, "end recip_lattice\n")
    @printf(io, "\n")
    return nothing
end

"""Write `kpoints` block."""
function _nnkp_write_block_kpoints(io::IO, kpoints)
    @printf(io, "begin kpoints\n")
    @printf(io, "%d\n", length(kpoints))
    for kpt in kpoints
        @printf(io, "%14.8f %14.8f %14.8f\n", kpt...)
    end
    @printf(io, "end kpoints\n")
    @printf(io, "\n")
    return nothing
end

"""Write `projections` block."""
function _nnkp_write_block_projections(io::IO, projections)
    @printf(io, "begin projections\n")
    @printf(io, "%d\n", length(projections))
    for p in projections
        @printf(io, "%12.7f %12.7f %12.7f %3d %3d %3d\n", p.center..., p.l, p.m, p.n)
        @printf(io, "%12.7f %12.7f %12.7f    ", p.zaxis...)
        @printf(io, "%12.7f %12.7f %12.7f    ", p.xaxis...)
        @printf(io, "%8.4f\n", p.α)
    end
    @printf(io, "end projections\n")
    @printf(io, "\n")
    return nothing
end

"""Write `auto_projections` block."""
function _nnkp_write_block_auto_projections(io::IO, auto_projections)
    @printf(io, "begin auto_projections\n")
    @printf(io, "%d\n", auto_projections)
    @printf(io, "%d\n", 0)  # magic number, useless
    @printf(io, "end auto_projections\n")
    @printf(io, "\n")
    return nothing
end

"""Write `nnkpts` block."""
function _nnkp_write_block_nnkpts(io::IO, kpb_k, kpb_G)
    n_kpts = length(kpb_k)
    n_bvecs = length(kpb_k[1])

    @printf(io, "begin nnkpts\n")
    @printf(io, "%d\n", n_bvecs)
    for ik in 1:n_kpts
        for ib in 1:n_bvecs
            @printf(io, "%6d %6d %6d %6d %6d\n", ik, kpb_k[ik][ib], kpb_G[ik][ib]...)
        end
    end
    @printf(io, "end nnkpts\n")
    @printf(io, "\n")
    return nothing
end

"""Write `exclude_bands` block."""
function _nnkp_write_block_exclude_bands(io::IO, exclude_bands)
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

function write_nnkp(io::IO, params::AbstractDict, ::W90InputToml; header = default_header())
    _nnkp_check_required_params(params)
    _check_dimensions_kpb(params["kpb_k"], params["kpb_G"])
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
        format::AbstractFileFormat;
        header = default_header(),
    )
    return open(filename, "w") do io
        write_nnkp(io, params, format; header)
    end
end

function write_nnkp(
        file::Union{IO, AbstractString}, params::AbstractDict; header = default_header()
    )
    format = w90input_format()
    return write_nnkp(file, params, format; header)
end
