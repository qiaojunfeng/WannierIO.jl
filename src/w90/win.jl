using Printf: @printf, @sprintf

export read_win, write_win

"""
    read_win(filename; fix_inputs=true)

Read the input file of `Wannier90`.

# Arguments
- `filename::AbstractString`: The name of the input file.

# Keyword Arguments
- `fix_inputs`: sanity check and fix the input parameters, e.g., set
    `num_bands = num_wann` if `num_bands` is not specified,
    convert `atoms_cart` always to `atoms_frac`, etc.
    See also [`_fix_win!`](@ref).
"""
function read_win(filename::AbstractString; fix_inputs::Bool=true)
    @info "Reading win file: $filename"
    io = open(filename)

    # The win file uses "num_wann", so I keep it as is, and not using "n_wann".
    keys_int = [
        :num_wann,
        :num_bands,
        :num_iter,
        :conv_window,
        :num_cg_steps,
        :wannier_plot_supercell,
    ]
    keys_int3 = [:mp_grid]
    keys_float = [
        :kmesh_tol,
        :dis_froz_min,
        :dis_froz_max,
        :dis_win_min,
        :dis_win_max,
        :fermi_energy,
        :conv_tol,
    ]
    keys_bool = [
        :use_ws_distance,
        :wannier_plot,
        :wvfn_formatted,
        :write_hr,
        :spn_formatted,
        :write_tb,
        :bands_plot,
    ]

    params = Dict{Symbol,Any}()

    read_line() = strip(readline(io))
    remove_comments(line::AbstractString) = begin
        i = findfirst(r"!|#", line)
        if i !== nothing
            line = strip(line[1:(i.start - 1)])
        end
        return line
    end
    # handle case insensitive win files (relic of Fortran)
    read_line_until_nonempty(; lower=true) = begin
        while !eof(io)
            line = read_line()
            lower && (line = lowercase(line))
            line = remove_comments(line)
            if !isempty(line)
                return line
            end
        end
        error("Unexpected end of file.")
    end
    parse_array(line::AbstractString; T=Float64) = map(x -> parse(T, x), split(line))
    read_array(f::IOStream) = parse_array(readline(f))

    while !eof(io)
        line = read_line_until_nonempty()

        # first handle special cases, e.g., blocks
        if occursin(r"begin\s+unit_cell_cart", line)
            unit_cell = zeros(Float64, 3, 3)
            unit = read_line_until_nonempty()
            if !startswith(unit, "b") && !startswith(unit, "a")
                line = unit
                unit = "ang"
            else
                line = read_line_until_nonempty()
            end
            for i in 1:3
                # in win file, each line is a lattice vector, here it is stored as column vec
                unit_cell[:, i] = parse_array(line)
                line = read_line_until_nonempty()
            end
            @assert occursin(r"end\s+unit_cell_cart", line)
            if startswith(unit, "b")
                # convert to angstrom
                unit_cell .*= Bohr
            end
            unit_cell = Mat3{Float64}(unit_cell)
            push!(params, :unit_cell_cart => unit_cell)
        elseif occursin(r"begin\s+atoms_(frac|cart)", line)
            iscart = occursin(r"cart", line)
            # do not lowercase due to atomic label
            line = read_line_until_nonempty(; lower=false)

            if iscart
                unit = lowercase(line)
                if !startswith(unit, "b") && !startswith(unit, "a")
                    unit = "ang"
                else
                    # do not lowercase due to atomic label
                    line = read_line_until_nonempty(; lower=false)
                end
            end

            # I need to read all lines and get n_atoms
            lines = Vector{String}()
            while !occursin(r"end\s+atoms_(frac|cart)", lowercase(line))
                push!(lines, line)
                line = read_line_until_nonempty(; lower=false)
            end
            n_atoms = length(lines)
            atoms_frac = zeros(Float64, 3, n_atoms)
            atom_labels = Vector{String}()
            for i in 1:n_atoms
                l = split(lines[i])
                push!(atom_labels, l[1])
                atoms_frac[:, i] = parse_float.(l[2:end])
            end

            if iscart
                if startswith(unit, "b")
                    # convert to angstrom
                    atoms_frac .*= Bohr
                end
                push!(params, :atoms_cart => atoms_frac)
            else
                push!(params, :atoms_frac => atoms_frac)
            end
            push!(params, :atom_labels => atom_labels)
        elseif occursin(r"begin\s+projections", line)
            projections = Vector{String}()
            line = read_line_until_nonempty(; lower=false)
            while !occursin(r"end\s+projections", lowercase(line))
                push!(projections, line)
                line = read_line_until_nonempty(; lower=false)
            end
            push!(params, :projections => projections)
        elseif occursin(r"begin\s+kpoints", line)
            line = read_line_until_nonempty()
            # I need to read all lines and get n_kpts
            lines = Vector{String}()
            while !occursin(r"end\s+kpoints", line)
                push!(lines, line)
                line = read_line_until_nonempty()
            end

            n_kpts = length(lines)
            kpoints = zeros(Float64, 3, n_kpts)
            for i in 1:n_kpts
                # There might be weight at 4th column, but we don't use it.
                kpoints[:, i] = parse_array(lines[i])[1:3]
            end
            push!(params, :kpoints => kpoints)
        elseif occursin(r"begin\skpoint_path", line)
            KV = Pair{Symbol,Vec3{Float64}}
            kpoint_path = Vector{Vector{KV}}()

            # allow uppercase
            line = read_line_until_nonempty(; lower=false)
            while !occursin(r"end\s+kpoint_path", lowercase(line))
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

                line = read_line_until_nonempty(; lower=false)
            end
            push!(params, :kpoint_path => kpoint_path)
        else
            # now treat remaining lines as key-value pairs
            line = strip(replace(line, "=" => " ", ":" => " ", "," => " "))
            key, value = split(line; limit=2)
            key = Symbol(key)
            if key in keys_int
                value = parse(Int, value)
            elseif key in keys_int3
                value = parse_array(value; T=Int)
            elseif key in keys_float
                value = parse_float(value)
            elseif key in keys_bool
                value = parse_bool(value)
            end
            push!(params, key => value)
        end
    end
    close(io)

    fix_inputs && _fix_win!(params)

    # convert to NamedTuple, easier to access its fields with dot notation,
    # e.g., params.num_wann
    params = NamedTuple(params)

    println("  num_wann  = ", params.num_wann)
    println("  num_bands = ", params.num_bands)
    @printf("  mp_grid   = %d %d %d\n", params.mp_grid...)
    println()

    return params
end

"""
    _fix_win!(params)

Sanity check and add missing input parameters from a `win` file.

See also [`read_win`](@ref).
"""
function _fix_win!(params::Dict)
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
        size(params[:kpoints]) != (3, n_kpts) && error("kpoints has wrong shape")
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
        atoms_frac = inv(params[:unit_cell_cart]) * atoms_cart
        push!(params, :atoms_frac => atoms_frac)
    end
    return nothing
end

"""
    write_win(filename::AbstractString; kwargs...)

Write input parameters into a Wannier90 `win` file.

The input parameters are keyword arguments, with key names same as that of
Wannier90. The only exception is `atoms_frac`, which contains both atom
labels and fractional coordinates; however, here it is split into two keywords:
`atom_labels` which is a vector of element names, and `atoms_frac` which is
a matrix.

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
    # both [:Si, :Si] and ["Si", "Si"] are allowed
    atom_labels = ["Si", "Si"],
    # atoms_frac is a matrix, its columns are the fractional coordinates
    atoms_frac=[
        0.0  0.25
        0.0  0.25
        0.0  0.25
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
        0.0  0.0  0.0  0.0  0.5  0.5  0.5  0.5
        0.0  0.0  0.5  0.5  0.0  0.0  0.5  0.5
        0.0  0.5  0.0  0.5  0.0  0.5  0.0  0.5
    ],
    # additional parameters can be passed as keyword arguments, e.g.,
    num_iter=500,
)
```
"""
function write_win(filename::AbstractString; kwargs...)
    required_keys = [
        :num_wann, :unit_cell_cart, :atom_labels, :atoms_frac, :mp_grid, :kpoints
    ]
    for k in required_keys
        haskey(kwargs, k) || error("Required parameter $k not found")
    end

    # convert immutable kwargs to Dict
    params = Dict(kwargs)
    num_wann = pop!(params, :num_wann)
    num_bands = pop!(params, :num_bands, nothing)
    unit_cell_cart = pop!(params, :unit_cell_cart)
    atom_labels = pop!(params, :atom_labels)
    atoms_frac = pop!(params, :atoms_frac)
    @assert length(atom_labels) == size(atoms_frac, 2)
    projections = pop!(params, :projections, nothing)
    kpoint_path = pop!(params, :kpoint_path, nothing)
    mp_grid = pop!(params, :mp_grid)
    kpoints = pop!(params, :kpoints)

    open(filename, "w") do fp
        println(fp, "# Generated by WannierIO.jl at $(now())\n")

        # First write most important parameters
        println(fp, "num_wann = $num_wann")
        !isnothing(num_bands) && println(fp, "num_bands = $num_bands")
        # a new line to separate from the rest
        println(fp)

        # Drop remaining parameters
        for (key, value) in pairs(params)
            ismissing(value) && continue
            # @printf fp "%-20s = %-30s\n" string(key) value
            @printf fp "%s = %s\n" string(key) value
        end
        println(fp)

        # Lattice vectors are in rows in Wannier90
        println(fp, "begin unit_cell_cart\nangstrom")
        # for floats, e.g., unit_cell, kpoints, I need to write enough digits
        # to avoid rounding errors in finding bvectors
        for vec in eachcol(unit_cell_cart)
            @printf fp "%14.8f  %14.8f  %14.8f\n" vec...
        end
        println(fp, "end unit_cell_cart\n")

        println(fp, "begin atoms_frac")
        for (element, position) in zip(atom_labels, eachcol(atoms_frac))
            if isa(element, Symbol)
                element = string(element)
            end
            @printf fp "%-3s  %14.8f  %14.8f  %14.8f\n" element position...
        end
        println(fp, "end atoms_frac\n")

        if !isnothing(projections)
            println(fp, "begin projections")
            for proj in projections
                println(fp, proj)
            end
            println(fp, "end projections\n")
        end

        if !isnothing(kpoint_path)
            println(fp, "begin kpoint_path")
            for segment in kpoint_path
                line = ""
                for (label, kpt) in segment
                    seg = @sprintf "%-3s  %14.8f  %14.8f  %14.8f    " string(label) kpt...
                    line *= seg
                end
                line = rstrip(line)
                println(fp, line)
            end
            println(fp, "end kpoint_path\n")
        end

        @printf fp "mp_grid = %d  %d  %d\n\n" mp_grid...

        println(fp, "begin kpoints")
        for kpt in eachcol(kpoints)
            @printf fp "%14.8f  %14.8f  %14.8f\n" kpt...
        end
        println(fp, "end kpoints")
    end
    return nothing
end
