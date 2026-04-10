using Printf
using OrderedCollections
using TOML

export read_win, write_win

"""
    read_win(file; standardize=true)
    read_win(file, ::W90InputText; standardize=true)
    read_win(file, ::W90InputToml; standardize=true)

Read wannier90 input `win` file.

# Arguments
- `file`: The name of the input file, or an `IO`.

# Keyword Arguments
- `standardize`: sanity check and fix the input parameters, e.g., set
    `num_bands = num_wann` if `num_bands` is not specified,
    convert `atoms_cart` always to `atoms_frac`, etc.
    See also [`standardize_win!`](@ref).
"""
function read_win end

function read_win(io::IO, ::W90InputText; standardize::Bool = true)
    params = OrderedDict{String, Any}()

    while !eof(io)
        line = nextline(io)
        isempty(line) && break

        isblock, content = _win_check_line(line)
        if isblock
            push!(params, _win_parse_block(io, content))
        else
            push!(params, _win_parse_keyval(content))
        end
    end

    _win_convert_keyval_types!(params)
    standardize && standardize_win!(params)

    return params
end

function read_win(io::IO, ::W90InputToml; standardize::Bool = true)
    win = read_toml(io)
    standardize && standardize_win!(win)
    return win
end

function read_win(
        filename::AbstractString, format::AbstractFileFormat; standardize::Bool = true
    )
    return open(filename) do io
        read_win(io, format; standardize)
    end
end

function read_win(file::Union{IO, AbstractString}; standardize::Bool = true)
    format = detect_w90input_format(file)
    win = read_win(file, format; standardize)
    return win
end

"""Read the next non-empty line from a win file block, with optional case control.

Win files are case-insensitive, but some blocks (e.g., atoms_frac) preserve
atomic labels. This function supports both via the `:lower` kwarg.
"""
function _win_block_nextline(io::IO, block_name; kwargs...)
    # win file is case-insensitive, I will lowercase lines by default, but
    # for some blocks, e.g., atoms_frac, I need to keep the original case for
    # atomic labels. So I allow kwargs to control this behavior.
    lower = get(kwargs, :lower, true)
    line = nextline(io; kwargs..., lower)
    if isempty(line)
        error("Error parsing block `$block_name`: unexpected end of file")
    end
    return line
end

"""Check if a line marks the end of a named block in a win file."""
@inline function _win_block_isend(
        line::AbstractString, block_name::AbstractString; lower::Bool = true
    )
    block_name = lower ? lowercase(block_name) : block_name
    line = lower ? lowercase(line) : line
    return occursin(r"^end\s+" * block_name, line)
end

"""Assert that a line marks the end of a named block, or raise an error."""
@inline function _win_block_mustend(
        line::AbstractString, block_name::AbstractString; lower::Bool = true
    )
    isend = _win_block_isend(line, block_name; lower)
    return isend || error("Error parsing $block_name: `end $block_name` not found")
end

"""Parse a `unit_cell_cart` block from a win file.

Returns lattice vectors as a 3×3 matrix (columns are lattice vectors in angstrom).
"""
function _win_parse_block_unit_cell_cart(io::IO)
    block_name = "unit_cell_cart"
    unit_cell = zeros(Float64, 3, 3)

    unit = _win_block_nextline(io, block_name)
    if !startswith(unit, "b") && !startswith(unit, "a")
        line = unit
        unit = "ang"
    else
        line = _win_block_nextline(io, block_name)
    end

    for i in 1:3
        # in win file, each line is a lattice vector, here it is stored as column vec
        unit_cell[:, i] = parse_vector(line)
        line = _win_block_nextline(io, block_name)
    end
    _win_block_mustend(line, block_name)

    if startswith(unit, "b")
        # convert to angstrom
        unit_cell .*= Bohr
    end
    return mat3(unit_cell)
end

"""Parse a `atoms_frac` or `atoms_cart` block from a win file.

Returns a vector of (atom_label => fractional/Cartesian coordinates) pairs.
"""
function _win_parse_block_atoms(io::IO, block_name::AbstractString)
    block_name == "atoms_frac" ||
        block_name == "atoms_cart" ||
        error("Invalid atoms block: $block_name")
    iscart = block_name == "atoms_cart"
    # do not lowercase due to atomic label
    line = _win_block_nextline(io, block_name; lower = false)

    if iscart
        unit = lowercase(line)
        if !startswith(unit, "b") && !startswith(unit, "a")
            unit = "ang"
        else
            # do not lowercase due to atomic label
            line = _win_block_nextline(io, block_name; lower = false)
        end
    end

    lines = String[]
    while !_win_block_isend(line, block_name)
        push!(lines, line)
        line = _win_block_nextline(io, block_name; lower = false)
    end

    atoms = StringVec3{Float64}[]
    for l in lines
        s = split(l)
        atom = s[1]
        frac = parse_vector(s[2:end])
        push!(atoms, stringvec3(atom, frac))
    end

    if iscart
        if startswith(unit, "b")
            # convert to angstrom
            atoms = map(atoms) do (atom, pos)
                stringvec3(atom, pos .* Bohr)
            end
        end
        return "atoms_cart" => atoms
    end
    return "atoms_frac" => atoms
end

"""Parse a `kpoints` block from a win file.

Returns a vector of k-point coordinates.
"""
function _win_parse_block_kpoints(io::IO)
    block_name = "kpoints"
    lines = String[]
    line = _win_block_nextline(io, block_name)
    while !_win_block_isend(line, block_name)
        push!(lines, line)
        line = _win_block_nextline(io, block_name)
    end

    kpoints = zeros(Vec3{Float64}, length(lines))
    for i in eachindex(lines)
        # There might be weight at 4th column, but we don't use it.
        kpoints[i] = vec3(parse_vector(lines[i])[1:3])
    end
    return kpoints
end

"""Parse a `kpoint_path` block from a win file.

Returns a vector of segments, each containing (start_label => start_kpt, end_label => end_kpt) pairs.
"""
function _win_parse_block_kpoint_path(io::IO)
    block_name = "kpoint_path"
    kpoint_path = Vector{Vector{StringVec3{Float64}}}()

    # allow uppercase
    line = _win_block_nextline(io, block_name; lower = false)
    while !_win_block_isend(line, block_name)
        l = split(line)
        length(l) == 8 || error("Invalid kpoint_path line: $line")
        # start kpoint
        start_label = l[1]
        start_kpt = vec3(parse_vector(l[2:4]))
        # end kpoint
        end_label = l[5]
        end_kpt = vec3(parse_vector(l[6:8]))
        push!(kpoint_path, [start_label => start_kpt, end_label => end_kpt])

        line = _win_block_nextline(io, block_name; lower = false)
    end
    return kpoint_path
end

"""Parse a `explicit_kpath` block from a win file.

Returns a vector of k-point coordinates along the path.
"""
function _win_parse_block_explicit_kpath(io::IO)
    block_name = "explicit_kpath"
    explicit_kpath = Vec3{Float64}[]

    line = _win_block_nextline(io, block_name)
    while !_win_block_isend(line, block_name)
        l = split(line)
        length(l) >= 3 || error("Invalid $block_name line: $line")
        push!(explicit_kpath, vec3(parse_vector(l[1:3])))
        line = _win_block_nextline(io, block_name)
    end
    return explicit_kpath
end

"""Parse a `explicit_kpath_labels` block from a win file.

Returns a vector of (label => k-point) pairs for high-symmetry points.
"""
function _win_parse_block_explicit_kpath_labels(io::IO)
    block_name = "explicit_kpath_labels"
    explicit_kpath_labels = StringVec3{Float64}[]

    # allow uppercase
    line = _win_block_nextline(io, block_name; lower = false)
    while !_win_block_isend(line, block_name)
        l = split(line)
        length(l) == 4 || error("Invalid $block_name line: $line")
        label = l[1]
        kpt = vec3(parse_vector(l[2:4]))
        push!(explicit_kpath_labels, label => kpt)
        line = _win_block_nextline(io, block_name; lower = false)
    end
    return explicit_kpath_labels
end

"""Parse a generic block as a vector of strings (for unknown block types)."""
function _win_parse_block_string(io::IO, block_name::AbstractString)
    block_content = String[]
    # allow uppercase
    line = _win_block_nextline(io, block_name; lower = false)
    while !_win_block_isend(line, block_name)
        push!(block_content, line)
        line = _win_block_nextline(io, block_name; lower = false)
    end
    return block_content
end

"""Dispatch to the appropriate block parser based on the block name."""
function _win_parse_block(io::IO, block_name::AbstractString)
    if block_name == "unit_cell_cart"
        return "unit_cell_cart" => _win_parse_block_unit_cell_cart(io)
    elseif block_name == "atoms_frac" || block_name == "atoms_cart"
        return _win_parse_block_atoms(io, block_name)
    elseif block_name == "projections"
        return "projections" => _win_parse_block_string(io, "projections")
    elseif block_name == "kpoints"
        return "kpoints" => _win_parse_block_kpoints(io)
    elseif block_name == "kpoint_path"
        return "kpoint_path" => _win_parse_block_kpoint_path(io)
    elseif block_name == "explicit_kpath"
        return "explicit_kpath" => _win_parse_block_explicit_kpath(io)
    elseif block_name == "explicit_kpath_labels"
        return "explicit_kpath_labels" => _win_parse_block_explicit_kpath_labels(io)
    end
    # Treat unknown blocks as Vector{String}.
    return block_name => _win_parse_block_string(io, block_name)
end

"""Return a dictionary mapping parameter names to their value types for parsing.

Used to determine how to parse key-value pairs from win files.
"""
function _win_keyval_types()
    key_types = Dict{String, Symbol}()

    for key in [
            "num_wann",
            "num_bands",
            "num_iter",
            "dis_num_iter",
            "dis_conv_window",
            "conv_window",
            "num_cg_steps",
            "num_print_cycles",
            "iprint",
            "search_shells",
            "bands_num_points",
            "ws_search_size",
            "num_guide_cycles",
            "num_no_guide_iter",
        ]
        key_types[key] = :int
    end

    for key in ["mp_grid", "wannier_plot_supercell"]
        key_types[key] = :int3
    end

    for key in [
            "kmesh_tol",
            "conv_tol",
            "dis_froz_min",
            "dis_froz_max",
            "dis_win_min",
            "dis_win_max",
            "dis_mix_ratio",
            "dis_conv_tol",
            "fermi_energy",
            "fermi_energy_min",
            "fermi_energy_max",
            "fermi_energy_step",
            "ws_distance_tol",
        ]
        key_types[key] = :float
    end

    for key in [
            "use_ws_distance",
            "wannier_plot",
            "bands_plot",
            "wvfn_formatted",
            "spn_formatted",
            "write_hr",
            "write_tb",
            "write_xyz",
            "write_rmn",
            "guiding_centres",
            "gamma_only",
            "spinors",
            "postproc_setup",
            "auto_projections",
            "restart",
        ]
        key_types[key] = :bool
    end

    for key in ["exclude_bands", "select_projections"]
        key_types[key] = :indices
    end

    return key_types
end

"""Parse a key-value line, separating key and value by = or : delimiters."""
function _win_parse_keyval(line::AbstractString)
    normalized_line = strip(replace(line, '=' => ' ', ':' => ' '))
    parts = split(normalized_line; limit = 2)
    length(parts) == 2 || error("Invalid key-value line: $line")

    key, value = parts
    return key => strip(value)  # remove leading whitespaces
end

"""Convert a string value to the appropriate type based on value_type symbol."""
function _win_convert_keyval_type(value::AbstractString, value_type::Symbol)
    if value_type == :int
        return parse(Int, value)
    elseif value_type == :int3
        parsed = parse_vector(strip(replace(value, "," => " ")), Int)
        return length(parsed) == 1 ? parsed[1] : parsed
    elseif value_type == :float
        return parse_float(value)
    elseif value_type == :bool
        return parse_bool(value)
    elseif value_type == :indices
        return parse_indices(value)
    end
    return value
end

"""Convert all key-value parameter strings to their appropriate types, in-place."""
function _win_convert_keyval_types!(params::AbstractDict)
    key_types = _win_keyval_types()
    for (key, value) in pairs(params)
        value_type = get(key_types, key, nothing)
        isnothing(value_type) && continue
        isa(value, AbstractString) || continue
        params[key] = _win_convert_keyval_type(value, value_type)
    end
    return params
end

"""Check if a line starts a block (begins with `'begin '`) or contains a key-value pair.

Returns `(isblock, content)` where `content` is the block name or key-value line.
"""
function _win_check_line(line::AbstractString)
    isblock = false
    content = line
    if startswith(line, "begin ")
        block_name = strip(line[7:end])
        isempty(block_name) && error("Invalid block line: $line")
        isblock = true
        # Force to String to be type stable
        content = string(block_name)
    end
    return isblock, content
end

"""
    $(SIGNATURES)

Sanity check and add missing input parameters from a `win` file.

See also [`read_win`](@ref).
"""
function standardize_win!(params::AbstractDict)
    !haskey(params, "num_wann") && error("num_wann not found")
    params["num_wann"] > 0 || throw(ArgumentError("num_wann must be positive"))

    # add num_bands if not found
    !haskey(params, "num_bands") && push!(params, "num_bands" => params["num_wann"])
    params["num_bands"] > 0 || throw(ArgumentError("num_bands must be positive"))

    if haskey(params, "mp_grid")
        length(params["mp_grid"]) != 3 && error("mp_grid has wrong length")
        any(i -> i <= 0, params["mp_grid"]) && error("mp_grid must be positive")
    else
        error("mp_grid not found")
    end

    if haskey(params, "kpoints")
        n_kpts = prod(params["mp_grid"])
        length(params["kpoints"]) != n_kpts && error("kpoints has wrong shape")
    else
        error("kpoints not found")
    end

    if haskey(params, "unit_cell_cart")
        any(x -> ismissing(x), params["unit_cell_cart"]) &&
            error("unit_cell_cart not found")
    else
        error("unit_cell_cart not found")
    end

    # if atoms_cart, convert to fractional
    if !haskey(params, "atoms_frac")
        !haskey(params, "atoms_cart") && error("both atoms_frac and atoms_cart are missing")
        atoms_cart = pop!(params, "atoms_cart")
        inv_cell = inv(params["unit_cell_cart"])
        atoms_frac = map(atoms_cart) do (atom, cart)
            stringvec3(atom, inv_cell * cart)
        end
        push!(params, "atoms_frac" => atoms_frac)
    end
    return nothing
end

"""
    write_win(file, params; header)
    write_win(file, params, ::W90InputText; header)
    write_win(file, params, ::W90InputToml; header)

Write input parameters into a wannier90 `win` file.

# Arguments
- `file`: The name of the output file, or an `IO`.
- `params`: a `Dict` (or `OrderedDict`) of parameters to be written into the `win` file

# Examples

```julia
using OrderedCollections, WannierIO

params = OrderedDict(
    "num_wann" => 4,
    "num_bands" => 4,
    # unit_cell_cart is a matrix, its columns are the lattice vectors in angstrom
    "unit_cell_cart" => [
        0.0      2.71527  2.71527
        2.71527  0.0      2.71527
        2.71527  2.71527  0.0
    ],
    # atoms_frac is a vector of pairs of atom_label and fractional coordinates
    "atoms_frac" => [
        "Si" => [0.0, 0.0, 0.0],
        "Si" => [0.25, 0.25, 0.25],
    ],
    # each element in projections will be written as a line in the win file
    "projections" => [
        "random",
    ],
    "kpoint_path" => [
        ["G" => [0.0, 0.0, 0.0], "X" => [0.5, 0.0, 0.5]],
        ["X" => [0.5, 0.0, 0.5], "U" => [0.625, 0.25, 0.625]],
    ],
    "mp_grid" => [2, 2, 2],
    # kpoints is a matrix, its columns are the fractional coordinates
    "kpoints" => [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.5],
        [0.0, 0.5, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.0],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
        [0.5, 0.5, 0.5],
    ],
    # additional parameters, e.g.,
    "num_iter" => 500,
)
write_win("silicon.win", params)
```
"""
function write_win end

function write_win(io::IO, params::AbstractDict, ::W90InputText; header = default_header())
    _win_check_required_params(params)

    # Copy params to an OrderedDict, to avoid modifying the input `params`.
    # Here we use OrderedDict to keep the order if the input is a OrderedDict.
    # Most likely the important parameters are at the beginning upon user input.
    params = OrderedDict(pairs(params))
    # These are blocks
    unit_cell_cart = pop!(params, "unit_cell_cart")
    atoms_frac = pop!(params, "atoms_frac", nothing)
    atoms_cart = pop!(params, "atoms_cart", nothing)
    projections = pop!(params, "projections", nothing)
    kpoint_path = pop!(params, "kpoint_path", nothing)
    explicit_kpath = pop!(params, "explicit_kpath", nothing)
    explicit_kpath_labels = pop!(params, "explicit_kpath_labels", nothing)
    kpoints = pop!(params, "kpoints", nothing)

    # an additional new line to separate from the rest
    isnothing(header) || _win_write_comment(io, header * "\n")
    _win_write_keyvals(io, params)
    # a new line to separate from the rest
    println(io)

    _win_write_block_unit_cell_cart(io, unit_cell_cart)
    isnothing(atoms_frac) || _win_write_block_atoms_frac(io, atoms_frac)
    isnothing(atoms_cart) || _win_write_block_atoms_cart(io, atoms_cart)
    isnothing(projections) || _win_write_block_projections(io, projections)
    isnothing(kpoint_path) || _win_write_block_kpoint_path(io, kpoint_path)
    isnothing(explicit_kpath_labels) ||
        _win_write_block_explicit_kpath_labels(io, explicit_kpath_labels)
    isnothing(explicit_kpath) || _win_write_block_explicit_kpath(io, explicit_kpath)
    isnothing(kpoints) || _win_write_block_kpoints(io, kpoints)

    return nothing
end

function write_win(io::IO, params::AbstractDict, ::W90InputToml; header = default_header())
    _win_check_required_params(params)
    isnothing(header) || println(io, header, "\n")
    write_toml(io, params)
    return nothing
end

function write_win(
        filename::AbstractString,
        params::AbstractDict,
        format::AbstractFileFormat;
        header = default_header(),
    )
    return open(filename, "w") do io
        write_win(io, params, format; header)
    end
end

function write_win(
        file::Union{IO, AbstractString}, params::AbstractDict; header = default_header()
    )
    return write_win(file, params, W90InputText(); header)
end

"""
    $(SIGNATURES)
"""
@inline function _win_check_required_params(kwargs)
    required_keys = ["num_wann", "unit_cell_cart", "mp_grid"]
    for k in required_keys
        haskey(kwargs, k) || throw(ArgumentError("Required parameter $k not found"))
    end
    atoms_cart_frac = haskey.(Ref(kwargs), ["atoms_cart", "atoms_frac"])
    return if all(atoms_cart_frac)
        error("Both atoms_cart and atoms_frac are found")
    elseif !any(atoms_cart_frac)
        error("Both atoms_frac and atoms_cart are missing")
    end
end

"""Write a comment line to the output, prefixing with # if needed."""
function _win_write_comment(io::IO, comment)
    return if !isnothing(comment)
        startswith(lstrip(comment), "#") || (comment = "# $comment")
        println(io, comment)
    end
end

"""Format a key-value pair for output, handling special types like int3 and indices."""
function _win_format_keyval(key::AbstractString, value, value_type::Union{Symbol, Nothing})
    if value_type == :int3
        # Unpack 3-element vector/tuple into separate integers
        return join(value, "  ")
    elseif value_type == :indices
        # Use format_indices for index arrays
        return isempty(value) ? "" : format_indices(value)
    else
        # Default: convert to string
        return string(value)
    end
end

"""Write all key-value parameters to the output stream."""
function _win_write_keyvals(io::IO, params::AbstractDict)
    key_types = _win_keyval_types()
    for (key, value) in pairs(params)
        ismissing(value) && continue
        # Skip empty indices
        if haskey(key_types, key) && key_types[key] == :indices && isempty(value)
            continue
        end

        value_type = get(key_types, key, nothing)
        formatted_value = _win_format_keyval(key, value, value_type)
        @printf io "%s = %s\n" string(key) formatted_value
    end
    return
end

"""Write a `unit_cell_cart` block to the output."""
function _win_write_block_unit_cell_cart(io::IO, unit_cell_cart)
    # Lattice vectors are in rows in Wannier90
    println(io, "begin unit_cell_cart\nangstrom")
    # for floats, e.g., unit_cell, kpoints, I need to write enough digits
    # to avoid rounding errors in finding bvectors
    for vec in eachcol(unit_cell_cart)
        @printf io "%14.8f  %14.8f  %14.8f\n" vec...
    end
    return println(io, "end unit_cell_cart\n")
end

"""Write a `atoms_frac` block to the output."""
function _win_write_block_atoms_frac(io::IO, atoms_frac)
    println(io, "begin atoms_frac")
    # label width is dynamic, e.g. 3 to accommodate "Si1"
    n = maximum(p -> length(string(first(p))), atoms_frac)
    fmt = Printf.Format("%-$(n)s  %14.8f  %14.8f  %14.8f\n")
    for (element, position) in atoms_frac
        Printf.format(io, fmt, string(element), position...)
    end
    return println(io, "end atoms_frac\n")
end

"""Write a `atoms_cart` block to the output."""
function _win_write_block_atoms_cart(io::IO, atoms_cart)
    println(io, "begin atoms_cart")
    # unit is angstrom
    println(io, "angstrom")
    # label width is dynamic, e.g. 3 to accommodate "Si1"
    n = maximum(p -> length(string(first(p))), atoms_cart)
    fmt = Printf.Format("%-$(n)s  %14.8f  %14.8f  %14.8f\n")
    for (element, position) in atoms_cart
        Printf.format(io, fmt, string(element), position...)
    end
    return println(io, "end atoms_cart\n")
end

"""Write a `projections` block to the output."""
function _win_write_block_projections(io::IO, projections)
    println(io, "begin projections")
    for proj in projections
        println(io, proj)
    end
    return println(io, "end projections\n")
end

"""Write a `kpoint_path` block to the output."""
function _win_write_block_kpoint_path(io::IO, kpoint_path)
    println(io, "begin kpoint_path")
    # label width is dynamic, e.g. 5 to accommodate "Gamma"
    n = maximum(maximum(p -> length(string(first(p))), seg) for seg in kpoint_path)
    fmt = Printf.Format("%-$(n)s  %14.8f  %14.8f  %14.8f    ")
    for segment in kpoint_path
        line = ""
        for (label, kpt) in segment
            seg = Printf.format(fmt, string(label), kpt...)
            line *= seg
        end
        line = rstrip(line)
        println(io, line)
    end
    return println(io, "end kpoint_path\n")
end

"""Write a `explicit_kpath_labels` block to the output."""
function _win_write_block_explicit_kpath_labels(io::IO, explicit_kpath_labels)
    println(io, "begin explicit_kpath_labels")
    # label width is dynamic, e.g. 5 to accommodate "Gamma"
    n = maximum(p -> length(string(first(p))), explicit_kpath_labels)
    fmt = Printf.Format("%-$(n)s  %14.8f  %14.8f  %14.8f\n")
    for (label, kpt) in explicit_kpath_labels
        Printf.format(io, fmt, string(label), kpt...)
    end
    return println(io, "end explicit_kpath_labels\n")
end

"""Write a `explicit_kpath` block to the output."""
function _win_write_block_explicit_kpath(io::IO, explicit_kpath)
    println(io, "begin explicit_kpath")
    for kpt in explicit_kpath
        @printf io "%14.8f  %14.8f  %14.8f\n" kpt...
    end
    return println(io, "end explicit_kpath\n")
end

"""Write a `kpoints` block to the output."""
function _win_write_block_kpoints(io::IO, kpoints)
    println(io, "begin kpoints")
    for kpt in kpoints
        @printf io "%14.8f  %14.8f  %14.8f\n" kpt...
    end
    return println(io, "end kpoints\n")
end
