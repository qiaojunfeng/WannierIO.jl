using TOML

export read_win, write_win

"""
    read_win(filename; fix_inputs=true)
    read_win(filename, ::Wannier90Text; fix_inputs=true)
    read_win(filename, ::Wannier90Toml; fix_inputs=true)

Read wannier90 input `win` file.

# Arguments
- `filename`: The name of the input file.

# Keyword Arguments
- `fix_inputs`: sanity check and fix the input parameters, e.g., set
    `num_bands = num_wann` if `num_bands` is not specified,
    convert `atoms_cart` always to `atoms_frac`, etc.
    See also [`fix_win!`](@ref).
"""
function read_win end

function read_win(filename::AbstractString, ::Wannier90Text; fix_inputs::Bool=true)
    params = open(filename) do io
        # The win file uses "num_wann", so I keep it as is, and not using "n_wann".
        keys_int = [
            :num_wann,
            :num_bands,
            :num_iter,
            :dis_num_iter,
            :dis_conv_window,
            :conv_window,
            :num_cg_steps,
            :wannier_plot_supercell,
            :num_print_cycles,
            :iprint,
            :search_shells,
            :bands_num_points,
            :ws_search_size,
            :num_guide_cycles,
            :num_no_guide_iter,
        ]
        keys_int3 = [:mp_grid]
        keys_float = [
            :kmesh_tol,
            :conv_tol,
            :dis_froz_min,
            :dis_froz_max,
            :dis_win_min,
            :dis_win_max,
            :dis_mix_ratio,
            :dis_conv_tol,
            :fermi_energy,
            :fermi_energy_min,
            :fermi_energy_max,
            :fermi_energy_step,
            :ws_distance_tol,
        ]
        keys_bool = [
            :use_ws_distance,
            :wannier_plot,
            :bands_plot,
            :wvfn_formatted,
            :spn_formatted,
            :write_hr,
            :write_tb,
            :write_xyz,
            :write_rmn,
            :guiding_centres,
            :gamma_only,
            :spinors,
            :postproc_setup,
            :auto_projections,
            :restart,
        ]
        keys_indices = [:exclude_bands, :select_projections]

        params = Dict{Symbol,Any}()

        read_line() = strip(readline(io))
        function remove_comments(line::AbstractString)
            i = findfirst(r"!|#", line)
            if i !== nothing
                line = strip(line[1:(i.start - 1)])
            end
            return line
        end
        # handle case insensitive win files (relic of Fortran)
        function read_line_until_nonempty(; lower=true, block_name=nothing)
            while !eof(io)
                line = read_line()
                lower && (line = lowercase(line))
                line = remove_comments(line)
                if !isempty(line)
                    return line
                end
            end
            # end of line reached
            if block_name !== nothing
                # in the middle of a block, raise error
                error("end of file reached while parsing block `$block_name`")
            end
            # else, outside of blocks, should be fine reaching end of file
            return nothing
        end
        function parse_array(line::AbstractString; T=Float64)
            return map(x -> parse(T, x), split(line))
        end
        read_array(f::IOStream) = parse_array(readline(f))

        while !eof(io)
            line = read_line_until_nonempty()
            line === nothing && break

            # first handle special cases, e.g., blocks
            if occursin(r"^begin\s+unit_cell_cart", line)
                block_name = "unit_cell_cart"
                unit_cell = zeros(Float64, 3, 3)
                unit = read_line_until_nonempty(; block_name)
                if !startswith(unit, "b") && !startswith(unit, "a")
                    line = unit
                    unit = "ang"
                else
                    line = read_line_until_nonempty(; block_name)
                end
                for i in 1:3
                    # in win file, each line is a lattice vector, here it is stored as column vec
                    unit_cell[:, i] = parse_array(line)
                    line = read_line_until_nonempty(; block_name)
                end
                @assert occursin(r"^end\s+unit_cell_cart", line) "error parsing $block_name: `end $block_name` not found"
                if startswith(unit, "b")
                    # convert to angstrom
                    unit_cell .*= Bohr
                end
                unit_cell = Mat3{Float64}(unit_cell)
                push!(params, :unit_cell_cart => unit_cell)
            elseif occursin(r"^begin\s+atoms_(frac|cart)", line)
                iscart = occursin("cart", line)
                block_name = "atoms_$(iscart ? "cart" : "frac")"
                # do not lowercase due to atomic label
                line = read_line_until_nonempty(; lower=false, block_name)

                if iscart
                    unit = lowercase(line)
                    if !startswith(unit, "b") && !startswith(unit, "a")
                        unit = "ang"
                    else
                        # do not lowercase due to atomic label
                        line = read_line_until_nonempty(; lower=false, block_name)
                    end
                end

                # I need to read all lines and get n_atoms
                lines = Vector{String}()
                while !occursin(Regex("^end\\s+" * block_name), lowercase(line))
                    push!(lines, line)
                    line = read_line_until_nonempty(; lower=false, block_name)
                end
                n_atoms = length(lines)
                atoms_frac = SymbolVec3{Float64}[]
                for i in 1:n_atoms
                    l = split(lines[i])
                    symbol = Symbol(l[1])
                    frac = Vec3(parse_float.(l[2:end])...)
                    push!(atoms_frac, SymbolVec3(symbol, frac))
                end

                if iscart
                    if startswith(unit, "b")
                        # convert to angstrom
                        atoms_frac = map(atoms_frac) do (symbol, pos)
                            SymbolVec3(symbol, pos .* Bohr)
                        end
                    end
                    push!(params, :atoms_cart => atoms_frac)
                else
                    push!(params, :atoms_frac => atoms_frac)
                end
            elseif occursin(r"^begin\s+projections", line)
                block_name = "projections"
                projections = Vector{String}()
                line = read_line_until_nonempty(; lower=false, block_name)
                while !occursin(r"^end\s+projections", lowercase(line))
                    push!(projections, line)
                    line = read_line_until_nonempty(; lower=false, block_name)
                end
                push!(params, :projections => projections)
            elseif occursin(r"^begin\s+kpoints", line)
                block_name = "kpoints"
                line = read_line_until_nonempty(; block_name)
                # I need to read all lines and get n_kpts
                lines = Vector{String}()
                while !occursin(r"^end\s+kpoints", line)
                    push!(lines, line)
                    line = read_line_until_nonempty(; block_name)
                end

                n_kpts = length(lines)
                kpoints = zeros(Vec3{Float64}, n_kpts)
                for i in 1:n_kpts
                    # There might be weight at 4th column, but we don't use it.
                    kpoints[i] = Vec3(parse_array(lines[i])[1:3])
                end
                push!(params, :kpoints => kpoints)
            elseif occursin(r"^begin\s+kpoint_path", line)
                block_name = "kpoint_path"
                kpoint_path = Vector{Vector{SymbolVec3{Float64}}}()

                # allow uppercase
                line = read_line_until_nonempty(; lower=false, block_name)
                while !occursin(r"^end\s+kpoint_path", lowercase(line))
                    l = split(line)
                    length(l) == 8 || error("Invalid kpoint_path line: $line")
                    # start kpoint
                    start_label = Symbol(l[1])
                    start_kpt = Vec3{Float64}(parse.(Float64, l[2:4]))
                    # end kpoint
                    end_label = Symbol(l[5])
                    end_kpt = Vec3{Float64}(parse.(Float64, l[6:8]))
                    # push to kpath
                    push!(kpoint_path, [start_label => start_kpt, end_label => end_kpt])

                    line = read_line_until_nonempty(; lower=false, block_name)
                end
                push!(params, :kpoint_path => kpoint_path)
            elseif occursin(r"^begin\s+(.+)", line)
                # treat all remaining unknown blocks as Vector of String
                block_name = match(r"^begin\s+(.+)", line).captures[1]
                block_content = Vector{String}()
                # allow uppercase
                line = read_line_until_nonempty(; lower=false, block_name)
                while !occursin(Regex("^end\\s+" * block_name), lowercase(line))
                    push!(block_content, line)
                    line = read_line_until_nonempty(; lower=false, block_name)
                end
                push!(params, Symbol(block_name) => block_content)
            else
                # now treat remaining lines as key-value pairs
                line = strip(replace(line, "=" => " ", ":" => " "))
                key, value = split(line; limit=2)
                value = strip(value)  # remove leading whitespaces
                key = Symbol(key)
                if key in keys_int
                    value = parse(Int, value)
                elseif key in keys_int3
                    value = strip(replace(value, "," => " "))
                    value = parse_array(value; T=Int)
                elseif key in keys_float
                    value = parse_float(value)
                elseif key in keys_bool
                    value = parse_bool(value)
                elseif key in keys_indices
                    value = parse_indices(value)
                end
                push!(params, key => value)
            end
        end
        return params
    end

    fix_inputs && fix_win!(params)

    # convert to NamedTuple, easier to access its fields with dot notation,
    # e.g., params.num_wann
    params = NamedTuple(params)
    return params
end

function read_win(filename::AbstractString, ::Wannier90Toml; fix_inputs::Bool=true)
    win = TOML.parsefile(filename)

    # I store atoms_frac and kpoint_path as Vector of SymbolVec3.
    # However, TOML.print does not accept Pair (specifically, SymbolVec3),
    # instead I convert SymbolVec3 to Dict in _write_win_toml.
    # On reading I convert it back.
    function convert_SymbolVec3(d::Dict)
        # SymbolVec3 are converted to Dict of length 1 when writing
        if length(d) == 1
            k, v = only(d)
            isa(k, String) && isa(v, Vector{<:Real}) && return SymbolVec3(k, v)
        end

        # Need to do the conversion recursively
        for (k, v) in pairs(d)
            if isa(v, Dict) || isa(v, Vector)
                d[k] = convert_SymbolVec3(v)
            end
        end
        return d
    end
    function convert_SymbolVec3(v::Vector)
        if isa(v, Vector{<:Dict}) || isa(v, Vector{<:Vector})
            return map(convert_SymbolVec3, v)
        end
        return v
    end
    win = convert_SymbolVec3(win)

    # Convert keys to Symbol
    win = Dict(Symbol(k) => v for (k, v) in pairs(win))

    fix_inputs && fix_win!(win)

    # Vector{Vector} -> Mat3
    if haskey(win, :unit_cell_cart)
        win[:unit_cell_cart] = Mat3(win[:unit_cell_cart])
    end

    # Vector{Vector} -> Vector{Vec3}
    if haskey(win, :kpoints)
        win[:kpoints] = [Vec3(k) for k in win[:kpoints]]
    end

    return NamedTuple(win)
end

function read_win(filename::AbstractString; fix_inputs=true)
    format = Wannier90Text()
    try
        TOML.parsefile(filename)
    catch err
        err isa TOML.ParserError || rethrow()
    else
        format = Wannier90Toml()
    end
    win = read_win(filename, format; fix_inputs)

    num_wann = win[:num_wann]
    num_bands = nothing
    if :num_bands in keys(win)
        num_bands = win[:num_bands]
    end
    # I need to convert to tuple so that @info does not output its type
    mp_grid = Tuple(win[:mp_grid])
    @info "Reading win file" filename num_wann num_bands mp_grid

    return win
end

"""
    $(SIGNATURES)

Sanity check and add missing input parameters from a `win` file.

See also [`read_win`](@ref).
"""
function fix_win!(params::Dict)
    !haskey(params, :num_wann) && error("num_wann not found")
    params[:num_wann] > 0 || error("num_wann must be positive")

    # add num_bands if not found
    !haskey(params, :num_bands) && push!(params, :num_bands => params[:num_wann])
    params[:num_bands] > 0 || error("num_bands must be positive")

    if haskey(params, :mp_grid)
        length(params[:mp_grid]) != 3 && error("mp_grid has wrong length")
        any(i -> i <= 0, params[:mp_grid]) && error("mp_grid must be positive")
    else
        error("mp_grid not found")
    end

    if haskey(params, :kpoints)
        n_kpts = prod(params[:mp_grid])
        length(params[:kpoints]) != n_kpts && error("kpoints has wrong shape")
    else
        error("kpoints not found")
    end

    if haskey(params, :unit_cell_cart)
        any(x -> ismissing(x), params[:unit_cell_cart]) && error("unit_cell_cart not found")
    else
        error("unit_cell_cart not found")
    end

    # if atoms_cart, convert to fractional
    if !haskey(params, :atoms_frac)
        !haskey(params, :atoms_cart) && error("both atoms_frac and atoms_cart are missing")
        atoms_cart = pop!(params, :atoms_cart)
        inv_cell = inv(params[:unit_cell_cart])
        atoms_frac = map(atoms_cart) do (symbol, cart)
            SymbolVec3{Float64}(symbol, inv_cell * cart)
        end
        push!(params, :atoms_frac => atoms_frac)
    end
    return nothing
end

"""
    write_win(filename; toml=false, kwargs...)
    write_win(filename, ::Wannier90Text; kwargs...)
    write_win(filename, ::Wannier90Toml; kwargs...)

Write input parameters into a wannier90 `win` file.

The input parameters are keyword arguments, with key names same as that of wannier90.

# Examples

```julia
using WannierIO

write_win(
    "silicon.win";
    num_wann=4,
    num_bands=4,
    # unit_cell_cart is a matrix, its columns are the lattice vectors in angstrom
    unit_cell_cart=[
        0.0      2.71527  2.71527
        2.71527  0.0      2.71527
        2.71527  2.71527  0.0
    ],
    # atoms_frac is a vector of pairs of atom_label and fractional coordinates
    atoms_frac=[
        :Si => [0.0, 0.0, 0.0],
        :Si => [0.25, 0.25, 0.25],
        # both `:Si` and `"Si"` are allowed
        # "Si" => [0.25, 0.25, 0.25],
    ],
    # each element in projections will be written as a line in the win file
    projections=[
        "random",
    ]
    kpoint_path=[
        [:G => [0.0, 0.0, 0.0], :X => [0.5, 0.0, 0.5]],
        [:X => [0.5, 0.0, 0.5], :U => [0.625, 0.25, 0.625]],
    ],
    mp_grid=[2, 2, 2],
    # kpoints is a matrix, its columns are the fractional coordinates
    kpoints=[
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.5],
        [0.0, 0.5, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.0],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
        [0.5, 0.5, 0.5],
    ],
    # additional parameters can be passed as keyword arguments, e.g.,
    num_iter=500,
)
```
"""
function write_win end

"""
    $(SIGNATURES)
"""
@inline function _check_win_required_params(kwargs)
    required_keys = [:num_wann, :unit_cell_cart, :atoms_frac, :mp_grid, :kpoints]
    for k in required_keys
        @assert haskey(kwargs, k) "Required parameter $k not found"
    end
end

function write_win(
    filename::AbstractString, ::Wannier90Text; header=default_header(), kwargs...
)
    _check_win_required_params(kwargs)

    # convert immutable kwargs to Dict
    params = Dict(kwargs)
    num_wann = pop!(params, :num_wann)
    num_bands = pop!(params, :num_bands, nothing)
    unit_cell_cart = pop!(params, :unit_cell_cart)
    atoms_frac = pop!(params, :atoms_frac)
    projections = pop!(params, :projections, nothing)
    kpoint_path = pop!(params, :kpoint_path, nothing)
    mp_grid = pop!(params, :mp_grid)
    kpoints = pop!(params, :kpoints)
    # if !startswith(lstrip(header), "#")
    #     header = "# $header"
    # end

    open(filename, "w") do io
        println(io, "$(header)\n")

        # First write most important parameters
        println(io, "num_wann = $num_wann")
        !isnothing(num_bands) && println(io, "num_bands = $num_bands")
        # a new line to separate from the rest
        println(io)

        # Drop remaining parameters
        for (key, value) in pairs(params)
            ismissing(value) && continue
            # @printf fp "%-20s = %-30s\n" string(key) value
            @printf io "%s = %s\n" string(key) value
        end
        println(io)

        # Lattice vectors are in rows in Wannier90
        println(io, "begin unit_cell_cart\nangstrom")
        # for floats, e.g., unit_cell, kpoints, I need to write enough digits
        # to avoid rounding errors in finding bvectors
        for vec in eachcol(unit_cell_cart)
            @printf io "%14.8f  %14.8f  %14.8f\n" vec...
        end
        println(io, "end unit_cell_cart\n")

        println(io, "begin atoms_frac")
        for (element, position) in atoms_frac
            if isa(element, Symbol)
                element = string(element)
            end
            @printf io "%-3s  %14.8f  %14.8f  %14.8f\n" element position...
        end
        println(io, "end atoms_frac\n")

        if !isnothing(projections)
            println(io, "begin projections")
            for proj in projections
                println(io, proj)
            end
            println(io, "end projections\n")
        end

        if !isnothing(kpoint_path)
            println(io, "begin kpoint_path")
            for segment in kpoint_path
                line = ""
                for (label, kpt) in segment
                    seg = @sprintf "%-3s  %14.8f  %14.8f  %14.8f    " string(label) kpt...
                    line *= seg
                end
                line = rstrip(line)
                println(io, line)
            end
            println(io, "end kpoint_path\n")
        end

        @printf io "mp_grid = %d  %d  %d\n\n" mp_grid...

        println(io, "begin kpoints")
        for kpt in kpoints
            @printf io "%14.8f  %14.8f  %14.8f\n" kpt...
        end
        println(io, "end kpoints")
    end
end

function write_win(
    filename::AbstractString, ::Wannier90Toml; header=default_header(), kwargs...
)
    _check_win_required_params(kwargs)

    open(filename, "w") do io
        println(io, header, "\n")
        write_toml(io; kwargs...)
    end
end

function write_win(filename::AbstractString; toml=false, kwargs...)
    if toml
        format = Wannier90Toml()
    else
        format = Wannier90Text()
    end

    num_wann = get(kwargs, :num_wann, nothing)
    num_bands = get(kwargs, :num_bands, nothing)
    mp_grid = get(kwargs, :mp_grid, nothing)
    # I need to convert to tuple so that @info does not output its type
    mp_grid !== nothing && (mp_grid = Tuple(mp_grid))
    @info "Writing win file" filename num_wann num_bands mp_grid

    return write_win(filename, format; kwargs...)
end
