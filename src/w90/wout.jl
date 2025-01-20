export read_wout

"""
    $(SIGNATURES)

Parse wannier90 `wout` file.

# Keyword Arguments
- `iterations`: if `true`, parse all the iterations of disentanglement and
    maximal localization. Default is `false`.

# Return
- `lattice`: each column is a lattice vector in Å
- `recip_lattice`: each column is a reciprocal lattice vector in Å⁻¹
- `atom_labels`: atomic symbols
- `atom_positions`: in fractional coordinates
- `kgrid`: kpoint grid used in Wannierization
- `centers`: center of each final WF, in Å
- `spreads`: spread of each final WF, in Å²
- `sum_centers`: sum of final WF centers, in Å
- `sum_spreads`: sum of final WF spreads, in Å²
- `ΩI`, `ΩD`, `ΩOD`, `Ωtotal`: final spread (components) in Å²
- `phase_factors`: optional, global (multiplicative) phase factor for obtaining
    real-valued (or close to real) MLWFs
- `im_re_ratios`: optional, maximum Im/Re ratio
- `iterations`: disentanglement and max localization convergence history,
    only parsed when kwarg `iterations=true`
"""
function read_wout(filename::AbstractString; iterations::Bool=false)
    return open(filename) do io
        mark_lattice = "Lattice Vectors ("
        mark_recip = "Reciprocal-Space Vectors ("
        mark_atom_start = "|   Site       Fractional Coordinate          Cartesian Coordinate"
        mark_atom_end = "*----------------------------------------------------------------------------*"
        mark_kgrid = "Grid size ="
        mark_finalstate_start = "Final State"
        mark_finalstate_end = "Sum of centres and spreads"
        mark_ΩI = "Omega I      ="
        mark_ΩD = "Omega D      ="
        mark_ΩOD = "Omega OD     ="
        mark_Ωtotal = "Omega Total  ="
        mark_phase = "Phase Factor ="
        mark_imre = "Maximum Im/Re Ratio ="
        # convergence history
        mark_dis_start = "Extraction of optimally-connected subspace"
        mark_dis_end = "Time to disentangle bands"
        mark_maxloc_start = "| Iter  Delta Spread     RMS Gradient      Spread (Ang^2)      Time  |<-- CONV"
        mark_maxloc_end = mark_finalstate_start
        # parsed results
        results = Dict{Symbol,Any}()
        iters = Dict{Symbol,Any}()

        srline() = strip(readline(io))

        while !eof(io)
            line = srline()
            if occursin(mark_lattice, line)
                lines = String[line]
                append!(lines, srline() for _ in 1:3)
                push!(results, :lattice => Mat3(_parse_wout_lattice(lines)))
                continue
            end
            if occursin(mark_recip, line)
                lines = String[line]
                append!(lines, srline() for _ in 1:3)
                push!(results, :recip_lattice => Mat3(_parse_wout_recip_lattice(lines)))
                continue
            end
            if occursin(mark_atom_start, line)
                lines = String[line]
                while line != mark_atom_end
                    line = srline()
                    push!(lines, line)
                end
                atom_labels, atom_positions = _parse_wout_atoms(lines)
                push!(
                    results, :atom_labels => atom_labels, :atom_positions => atom_positions
                )
                continue
            end
            if occursin(mark_kgrid, line)
                #  parse line `Grid size =  9 x  9 x  9      Total points =  729`
                line = split(line)[[4, 6, 8]]
                push!(results, :kgrid => parse.(Int, line))
                continue
            end
            if iterations
                if occursin(mark_dis_start, line)
                    lines = String[line]
                    while !occursin(mark_dis_end, line)
                        line = srline()
                        push!(lines, line)
                    end
                    push!(iters, :disentangle => _parse_wout_disentangle(lines))
                    continue
                end
                if occursin(mark_maxloc_start, line)
                    lines = String[line]
                    while !occursin(mark_maxloc_end, line)
                        line = srline()
                        push!(lines, line)
                    end
                    push!(iters, :wannierize => _parse_wout_wannierize(lines))
                    # Do not continue: we use the same "Final State" mark for
                    # checking both the end of max localization and the start of
                    # final spread. Therefore, we need to proceed with the next
                    # if block on parsing final state.
                end
            end
            if occursin(mark_finalstate_start, line)
                line = srline()
                lines = String[line]
                while !occursin(mark_finalstate_end, line)
                    line = srline()
                    push!(lines, line)
                end
                c, s, sc, ss = _parse_wout_wf_center_spread(lines)
                push!(
                    results,
                    :centers => c,
                    :spreads => s,
                    :sum_centers => sc,
                    :sum_spreads => ss,
                )
                continue
            end
            if occursin(mark_ΩI, line)
                push!(results, :ΩI => parse_float(split(line, mark_ΩI)[2]))
                line = srline()
                push!(results, :ΩD => parse_float(split(line, mark_ΩD)[2]))
                line = srline()
                push!(results, :ΩOD => parse_float(split(line, mark_ΩOD)[2]))
                line = srline()
                push!(results, :Ωtotal => parse_float(split(line, mark_Ωtotal)[2]))
                continue
            end
            if occursin(mark_phase, line)
                phases = ComplexF64[]
                while occursin(mark_phase, line)
                    s = split(line, "=")
                    v = parse(ComplexF64, s[end])
                    push!(phases, v)
                    line = srline()  # there is an empty line after phase block
                end
                push!(results, :phase_factors => phases)
                continue
            end
            if occursin(mark_imre, line)
                imre = Float64[]
                while occursin(mark_imre, line)
                    s = split(line, "=")
                    v = parse(Float64, s[end])
                    push!(imre, v)
                    line = srline()  # there is an empty line after im/re block
                end
                push!(results, :im_re_ratios => imre)
                continue
            end
        end
        iterations && push!(results, :iterations => NamedTuple(iters))
        return NamedTuple(results)
    end
end

"""
Parse block
```
Lattice Vectors (Ang)
a_1     0.000000   2.715265   2.715265
a_2     2.715265   0.000000   2.715265
a_3     2.715265   2.715265   0.000000
```
"""
function _parse_wout_lattice(lines)
    @assert length(lines) == 4 "wrong number of lines for lattice"
    ang_unit = occursin("Lattice Vectors (Ang)", lines[1])
    @assert ang_unit "wout unit is not Angstrom, not supported yet"
    lattice = zeros(Float64, 3, 3)

    for (i, line) in enumerate(lines[2:end])
        line = split(line)
        @assert line[1] == "a_$i"
        lattice[:, i] = parse_float.(line[2:end])
    end
    return lattice
end

"""
Parse block
```
Reciprocal-Space Vectors (Ang^-1)
b_1    -1.157011   1.157011   1.157011
b_2     1.157011  -1.157011   1.157011
b_3     1.157011   1.157011  -1.157011
```
"""
function _parse_wout_recip_lattice(lines)
    @assert length(lines) == 4 "wrong number of lines for reciprocal lattice"
    @assert occursin("Reciprocal-Space Vectors", lines[1])
    recip_lattice = zeros(Float64, 3, 3)

    for (i, line) in enumerate(lines[2:end])
        line = split(line)
        @assert line[1] == "b_$i"
        recip_lattice[:, i] = parse_float.(line[2:end])
    end
    return recip_lattice
end

"""
Parse block
```
|   Site       Fractional Coordinate          Cartesian Coordinate (Ang)     |
+----------------------------------------------------------------------------+
| Si   1   0.00000   0.00000   0.00000   |    0.00000   0.00000   0.00000    |
| Si   2   0.25000   0.25000   0.25000   |    1.35763   1.35763   1.35763    |
*----------------------------------------------------------------------------*
```
"""
function _parse_wout_atoms(lines)
    lines = lines[3:(end - 1)]
    n_atom = length(lines)
    atom_labels = Vector{String}()
    atom_positions = zeros(Vec3{Float64}, n_atom)
    for (i, line) in enumerate(lines)
        line = split(line)
        @assert line[1] == "|" line
        push!(atom_labels, line[2])
        # cartesian
        # atom_positions[i] = Vec3(parse_float.(line[8:10]))
        # fractional
        atom_positions[i] = Vec3(parse_float.(line[4:6]))
    end
    return atom_labels, atom_positions
end

"""
Parse block
```
                  Extraction of optimally-connected subspace
                  ------------------------------------------
+---------------------------------------------------------------------+<-- DIS
|  Iter     Omega_I(i-1)      Omega_I(i)      Delta (frac.)    Time   |<-- DIS
+---------------------------------------------------------------------+<-- DIS
      1      25.38943399      21.32896063       1.904E-01      0.00    <-- DIS
      2      21.53095611      20.16097533       6.795E-02      0.01    <-- DIS
      3      20.40788223      19.35260423       5.453E-02      0.01    <-- DIS
      4      19.53883989      18.75563591       4.176E-02      0.02    <-- DIS
...
    341      16.22884440      16.22884440      -1.883E-10      2.43    <-- DIS
    342      16.22884440      16.22884440      -1.799E-10      2.44    <-- DIS

            <<<      Delta < 2.000E-10  over  3 iterations     >>>
            <<< Disentanglement convergence criteria satisfied >>>

        Final Omega_I    16.22884440 (Ang^2)

+----------------------------------------------------------------------------+

Time to disentangle bands      2.546 (sec)
```
"""
function _parse_wout_disentangle(lines)
    mark_iter = "<-- DIS"
    iter = Int[]
    ΩI_previous = Float64[]
    ΩI_current = Float64[]
    ΔΩI = Float64[]

    i_start = 6
    for line in lines[i_start:end]
        occursin(mark_iter, line) || continue

        line = split(line)
        push!(iter, parse(Int, line[1]))
        push!(ΩI_previous, parse(Float64, line[2]))
        push!(ΩI_current, parse(Float64, line[3]))
        push!(ΔΩI, parse(Float64, line[4]))
    end
    return (; iter, ΩI_previous, ΩI_current, ΔΩI)
end

"""
Parse block
```
WF centre and spread    1  ( -0.000005,  0.000021,  0.000023 )     2.56218734
WF centre and spread    2  (  0.000013, -0.000054,  0.000016 )     3.19493515
WF centre and spread    3  ( -0.000005, -0.000054, -0.000055 )     3.19482997
WF centre and spread    4  (  0.000012,  0.000015, -0.000058 )     3.19526437
WF centre and spread    5  (  1.357637,  1.357611,  1.357610 )     2.56218214
WF centre and spread    6  (  1.357619,  1.357684,  1.357617 )     3.19532825
WF centre and spread    7  (  1.357638,  1.357687,  1.357686 )     3.19513205
WF centre and spread    8  (  1.357620,  1.357617,  1.357694 )     3.19460833
Sum of centres and spreads (  5.430528,  5.430529,  5.430534 )    24.29446759
```
"""
function _parse_wout_wf_center_spread(lines)
    mark_wf = "WF centre and spread"
    mark_sum = "Sum of centres and spreads"
    centers = Vec3{Float64}[]
    spreads = Float64[]
    sum_centers = nothing
    sum_spreads = nothing
    for line in lines
        if occursin(mark_wf, line)
            sline = split(line, r"[ ,()]"; keepempty=false)
            idx = parse(Int, sline[5])
            push!(centers, Vec3(parse_float.(sline[6:8])))
            push!(spreads, parse_float(sline[9]))
            (idx != length(centers)) && error("WF index mismatch at $line")
        end
        if occursin(mark_sum, line)
            sline = split(line, r"[ ,()]"; keepempty=false)
            sum_centers = Vec3(parse_float.(sline[6:8]))
            sum_spreads = parse_float(sline[9])
        end
    end
    return centers, spreads, sum_centers, sum_spreads
end

"""
Parse block
```
*------------------------------- WANNIERISE ---------------------------------*
+--------------------------------------------------------------------+<-- CONV
| Iter  Delta Spread     RMS Gradient      Spread (Ang^2)      Time  |<-- CONV
+--------------------------------------------------------------------+<-- CONV

------------------------------------------------------------------------------
Initial State
 WF centre and spread    1  ( -0.000005,  0.000021,  0.000023 )     2.56218734
 WF centre and spread    2  (  0.000013, -0.000054,  0.000016 )     3.19493515
 WF centre and spread    3  ( -0.000005, -0.000054, -0.000055 )     3.19482997
 WF centre and spread    4  (  0.000012,  0.000015, -0.000058 )     3.19526437
 WF centre and spread    5  (  1.357637,  1.357611,  1.357610 )     2.56218214
 WF centre and spread    6  (  1.357619,  1.357684,  1.357617 )     3.19532825
 WF centre and spread    7  (  1.357638,  1.357687,  1.357686 )     3.19513205
 WF centre and spread    8  (  1.357620,  1.357617,  1.357694 )     3.19460833
 Sum of centres and spreads (  5.430528,  5.430529,  5.430534 )    24.29446759

     0     0.243E+02     0.0000000000       24.2944680346       2.48  <-- CONV
       O_D=      0.2135529 O_OD=      7.8520707 O_TOT=     24.2944680 <-- SPRD
------------------------------------------------------------------------------
Cycle:      1
 WF centre and spread    1  ( -0.000005,  0.000020,  0.000022 )     2.46316318
 WF centre and spread    2  (  0.000014, -0.000057,  0.000015 )     3.19187236
 WF centre and spread    3  ( -0.000005, -0.000057, -0.000058 )     3.19179103
 WF centre and spread    4  (  0.000012,  0.000014, -0.000061 )     3.19222621
 WF centre and spread    5  (  1.357637,  1.357612,  1.357611 )     2.46315800
 WF centre and spread    6  (  1.357618,  1.357687,  1.357618 )     3.19226713
 WF centre and spread    7  (  1.357637,  1.357690,  1.357689 )     3.19209166
 WF centre and spread    8  (  1.357619,  1.357619,  1.357697 )     3.19156919
 Sum of centres and spreads (  5.430528,  5.430529,  5.430534 )    24.07813875

     1    -0.216E+00     0.2558717278       24.0781391952       2.49  <-- CONV
       O_D=      0.2218113 O_OD=      7.6274835 O_TOT=     24.0781392 <-- SPRD
Delta: O_D=  0.8258326E-02 O_OD= -0.2245872E+00 O_TOT= -0.2163288E+00 <-- DLTA
------------------------------------------------------------------------------
Cycle:      2
...
------------------------------------------------------------------------------
Cycle:     45
WF centre and spread    1  (  0.000001,  0.000006,  0.000006 )     1.95373328
WF centre and spread    2  (  0.000016, -0.000065,  0.000019 )     3.27910139
WF centre and spread    3  ( -0.000008, -0.000065, -0.000066 )     3.27921479
WF centre and spread    4  (  0.000014,  0.000016, -0.000070 )     3.27965818
WF centre and spread    5  (  1.357631,  1.357627,  1.357627 )     1.95372427
WF centre and spread    6  (  1.357616,  1.357695,  1.357615 )     3.27949625
WF centre and spread    7  (  1.357641,  1.357699,  1.357697 )     3.27951005
WF centre and spread    8  (  1.357617,  1.357616,  1.357707 )     3.27899768
Sum of centres and spreads (  5.430528,  5.430529,  5.430534 )    23.58343588

   45    -0.186E-09     0.0000077396       23.5834363262       2.88  <-- CONV
      O_D=      0.2612087 O_OD=      7.0933833 O_TOT=     23.5834363 <-- SPRD
Delta: O_D=  0.2557594E-06 O_OD= -0.2559458E-06 O_TOT= -0.1863789E-09 <-- DLTA
------------------------------------------------------------------------------

           <<<     Delta < 2.000E-10  over  3 iterations     >>>
           <<< Wannierisation convergence criteria satisfied >>>
```
"""
function _parse_wout_wannierize(lines)
    mark_init = "Initial State"
    mark_cycle = "Cycle:"
    mark_sprd = "<-- SPRD"

    iter = Int[]
    centers = Vector{Vector{Vec3{Float64}}}()
    spreads = Vector{Vector{Float64}}()
    sum_centers = Vector{Vec3{Float64}}()
    sum_spreads = Float64[]
    ΩD = Float64[]
    ΩOD = Float64[]
    Ωtotal = Float64[]

    """
    i_start: line number of "WF centre and spread    1"
    i_end: line number of "<-- SPRD"
    """
    function append_iter!(i_start, i_end)
        (i_end > i_start) || error("Invalid range $i_start:$i_end")
        c, s, sc, ss = _parse_wout_wf_center_spread(lines[i_start:(i_end - 3)])
        push!(centers, c)
        push!(spreads, s)
        push!(sum_centers, sc)
        push!(sum_spreads, ss)
        sline = split(lines[i_end])
        push!(ΩD, parse_float(sline[2]))
        push!(ΩOD, parse_float(sline[4]))
        push!(Ωtotal, parse_float(sline[6]))
        return nothing
    end

    # Initial State
    i_start = i_end = -1
    for (i, line) in enumerate(lines)
        if occursin(mark_init, line)
            i_start = i + 1
            continue
        end
        if occursin(mark_sprd, line)
            i_end = i
            break
        end
    end
    (i_start > 0) || error("No `$mark_init` found")
    (i_end > i_start) || error("No `$mark_sprd` found")
    push!(iter, 0)
    append_iter!(i_start, i_end)

    # Remaining cycles
    # for each iteration: from `Cycle` to horizontal line, inclusive
    n_lines_block = i_end - i_start + 4
    i = i_end + 1
    while i <= length(lines)
        line = lines[i]
        if occursin(mark_cycle, line)
            idx = parse(Int, split(line)[2])
            push!(iter, idx)
            i_start = i + 1
            i_end = i + n_lines_block - 3
            append_iter!(i_start, i_end)
            i += n_lines_block
            continue
        end
        i += 1
    end
    return (; iter, centers, spreads, sum_centers, sum_spreads, ΩD, ΩOD, Ωtotal)
end
