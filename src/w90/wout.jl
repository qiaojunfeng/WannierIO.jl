export read_wout

const WOUT_MARKS = (;
    lattice = "Lattice Vectors (",
    recip = "Reciprocal-Space Vectors (",
    atom_start = "|   Site       Fractional Coordinate          Cartesian Coordinate",
    atom_end = "*----------------------------------------------------------------------------*",
    kgrid = "Grid size =",
    finalstate_start = "Final State",
    finalstate_end = "Sum of centres and spreads",
    ΩI = "Omega I      =",
    ΩD = "Omega D      =",
    ΩOD = "Omega OD     =",
    Ωtotal = "Omega Total  =",
    phase = "Phase Factor =",
    imre = "Maximum Im/Re Ratio =",
    dis_start = "Extraction of optimally-connected subspace",
    dis_end = "Time to disentangle bands",
    dis_iter = "<-- DIS",
    maxloc_start = "| Iter  Delta Spread     RMS Gradient      Spread (Ang^2)      Time  |<-- CONV",
    maxloc_end = "Final State",
    maxloc_wf_c_s = "WF centre and spread",
    maxloc_sum_c_s = "Sum of centres and spreads",
)

"""
    $(SIGNATURES)

Parse wannier90 `wout` file.

# Keyword Arguments
- `iterations`: if `true`, parse all the iterations of disentanglement and
    maximal localization. Default is `false`.

# Return
- `lattice`: each column is a lattice vector in Å
- `recip_lattice`: each column is a reciprocal lattice vector in Å⁻¹
- `atom_labels`: atomic labels
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
function read_wout(io::IO; iterations::Bool = false)
    # parsed results
    results = OrderedDict{String, Any}()
    iters = OrderedDict{String, Any}()

    while !eof(io)
        line = readstrip(io)
        line_kind = _wout_line_kind(line; iterations)

        if line_kind == :lattice
            results["lattice"] = _wout_parse_lattice(io)
        elseif line_kind == :recip_lattice
            results["recip_lattice"] = _wout_parse_recip_lattice(io)
        elseif line_kind == :atoms
            atom_labels, atom_positions = _wout_parse_atoms(io)
            results["atom_labels"] = atom_labels
            results["atom_positions"] = atom_positions
        elseif line_kind == :kgrid
            results["kgrid"] = _wout_parse_kgrid(line)
        elseif line_kind == :disentangle
            iters["disentangle"] = _wout_parse_disentangle(io)
        elseif line_kind == :wannierize
            iters["wannierize"] = _wout_parse_wannierize(io)
            # _wout_parse_wannierize use "Final State" as the marker,
            # this line is already consumed, we need to parse final state
            # immediately after that.
            wf = _wout_parse_final_state(io)
            results["centers"] = wf.centers
            results["spreads"] = wf.spreads
            results["sum_centers"] = wf.sum_centers
            results["sum_spreads"] = wf.sum_spreads
        elseif line_kind == :final_state
            wf = _wout_parse_final_state(io)
            results["centers"] = wf.centers
            results["spreads"] = wf.spreads
            results["sum_centers"] = wf.sum_centers
            results["sum_spreads"] = wf.sum_spreads
        elseif line_kind == :omega
            omega = _wout_parse_Ω(io, line)
            results["ΩI"] = omega.ΩI
            results["ΩD"] = omega.ΩD
            results["ΩOD"] = omega.ΩOD
            results["Ωtotal"] = omega.Ωtotal
        elseif line_kind == :phase
            results["phase_factors"] = _wout_parse_phase_factor(io, line)
        elseif line_kind == :imre
            results["im_re_ratios"] = _wout_parse_im_re_ratio(io, line)
        end
    end
    iterations && push!(results, "iterations" => iters)
    return results
end

function read_wout(filename::AbstractString; iterations::Bool = false)
    return open(filename) do io
        read_wout(io; iterations)
    end
end

function _wout_line_kind(line::AbstractString; iterations::Bool = false)
    if occursin(WOUT_MARKS.lattice, line)
        return :lattice
    elseif occursin(WOUT_MARKS.recip, line)
        return :recip_lattice
    elseif occursin(WOUT_MARKS.atom_start, line)
        return :atoms
    elseif occursin(WOUT_MARKS.kgrid, line)
        return :kgrid
    elseif iterations && occursin(WOUT_MARKS.dis_start, line)
        return :disentangle
    elseif iterations && occursin(WOUT_MARKS.maxloc_start, line)
        return :wannierize
    elseif occursin(WOUT_MARKS.finalstate_start, line)
        return :final_state
    elseif occursin(WOUT_MARKS.ΩI, line)
        return :omega
    elseif occursin(WOUT_MARKS.phase, line)
        return :phase
    elseif occursin(WOUT_MARKS.imre, line)
        return :imre
    end
    return :none
end

"""
Parse line `Grid size =  9 x  9 x  9      Total points =  729`
"""
function _wout_parse_kgrid(line::AbstractString)
    sline = split(line)[[4, 6, 8]]
    return parse_vector(sline, Int)
end

"""
`line` is the 1st, remaining lines will be read from `io`.
```
         Spreads (Ang^2)       Omega I      =     3.956862958
        ================       Omega D      =     0.008030049
                               Omega OD     =     0.501987969
    Final Spread (Ang^2)       Omega Total  =     4.466880976
```
"""
function _wout_parse_Ω(io::IO, line::AbstractString)
    ΩI = parse_float(split(line, WOUT_MARKS.ΩI)[2])
    line = readstrip(io)
    ΩD = parse_float(split(line, WOUT_MARKS.ΩD)[2])
    line = readstrip(io)
    ΩOD = parse_float(split(line, WOUT_MARKS.ΩOD)[2])
    line = readstrip(io)
    Ωtotal = parse_float(split(line, WOUT_MARKS.Ωtotal)[2])
    return (; ΩI, ΩD, ΩOD, Ωtotal)
end

"""
Parse blocks like
```
      Wannier Function Num:    1       Phase Factor =    0.996157  +0.087588i
      Wannier Function Num:    2       Phase Factor =    0.996157  +0.087588i
      Wannier Function Num:    3       Phase Factor =    0.996157  +0.087588i
      Wannier Function Num:    4       Phase Factor =    0.998869  +0.047543i

```
or
```
      Wannier Function Num:    1       Maximum Im/Re Ratio =    4.566451
      Wannier Function Num:    2       Maximum Im/Re Ratio =    4.566481
      Wannier Function Num:    3       Maximum Im/Re Ratio =    4.566335
      Wannier Function Num:    4       Maximum Im/Re Ratio =    2.154381

```
The `line` contains the 1st, the remaining lines will be read from `io`.
The last line should be an empty line.
"""
function _wout_parse_repeated_equals(
        io::IO, line::AbstractString, marker::AbstractString, ::Type{T}
    ) where {T}
    values = T[]
    while occursin(marker, line)
        s = split(line, "=")
        push!(values, parse(T, s[end]))
        line = readstrip(io)  # there is an empty line after this block
    end
    return values
end

"""
Parse block
```
      Wannier Function Num:    1       Phase Factor =    0.996157  +0.087588i
      Wannier Function Num:    2       Phase Factor =    0.996157  +0.087588i
      Wannier Function Num:    3       Phase Factor =    0.996157  +0.087588i
      Wannier Function Num:    4       Phase Factor =    0.998869  +0.047543i

```
"""
function _wout_parse_phase_factor(io::IO, line::AbstractString)
    return _wout_parse_repeated_equals(io, line, WOUT_MARKS.phase, ComplexF64)
end

"""
Parse block
```
      Wannier Function Num:    1       Maximum Im/Re Ratio =    4.566451
      Wannier Function Num:    2       Maximum Im/Re Ratio =    4.566481
      Wannier Function Num:    3       Maximum Im/Re Ratio =    4.566335
      Wannier Function Num:    4       Maximum Im/Re Ratio =    2.154381

```
"""
function _wout_parse_im_re_ratio(io::IO, line::AbstractString)
    return _wout_parse_repeated_equals(io, line, WOUT_MARKS.imre, Float64)
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
function _wout_parse_lattice(io::IO)
    lattice = zeros(Float64, 3, 3)

    for i in 1:3
        parts = split(readstrip(io))
        parts[1] == "a_$i" || error("line does not start with a_$i")
        lattice[:, i] = parse_float.(parts[2:end])
    end
    return mat3(lattice)
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
function _wout_parse_recip_lattice(io::IO)
    recip_lattice = zeros(Float64, 3, 3)

    for i in 1:3
        parts = split(readstrip(io))
        parts[1] == "b_$i" || error("line does not start with b_$i")
        recip_lattice[:, i] = parse_float.(parts[2:end])
    end
    return mat3(recip_lattice)
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
function _wout_parse_atoms(io::IO)
    readstrip(io) # separator line

    atom_labels = String[]
    atom_positions = Vec3{Float64}[]
    while !eof(io)
        line = readstrip(io)
        line == WOUT_MARKS.atom_end && break
        parts = split(line)
        parts[1] == "|" || error("line does not start with |")
        push!(atom_labels, parts[2])
        push!(atom_positions, vec3(parse_float.(parts[4:6])))
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
function _wout_parse_disentangle(io::IO)
    iter = Int[]
    ΩI_previous = Float64[]
    ΩI_current = Float64[]
    ΔΩI = Float64[]

    while !eof(io)
        line = readstrip(io)
        occursin(WOUT_MARKS.dis_end, line) && break
        occursin(WOUT_MARKS.dis_iter, line) || continue

        parts = split(line)
        isnothing(tryparse(Int, parts[1])) && continue
        push!(iter, parse(Int, parts[1]))
        push!(ΩI_previous, parse(Float64, parts[2]))
        push!(ΩI_current, parse(Float64, parts[3]))
        push!(ΔΩI, parse(Float64, parts[4]))
    end
    return OrderedDict{String, Any}(
        "iter" => iter,
        "ΩI_previous" => ΩI_previous,
        "ΩI_current" => ΩI_current,
        "ΔΩI" => ΔΩI,
    )
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
function _wout_parse_wf_center_spread(io::IO)
    centers = Vec3{Float64}[]
    spreads = Float64[]
    sum_centers = nothing
    sum_spreads = nothing
    line = readstrip(io)

    while true
        if occursin(WOUT_MARKS.maxloc_wf_c_s, line)
            sline = split(line, r"[ ,()]"; keepempty = false)
            idx = parse(Int, sline[5])
            push!(centers, Vec3(parse_float.(sline[6:8])))
            push!(spreads, parse_float(sline[9]))
            (idx != length(centers)) && error("WF index mismatch at $line")
        end
        if occursin(WOUT_MARKS.maxloc_sum_c_s, line)
            sline = split(line, r"[ ,()]"; keepempty = false)
            sum_centers = Vec3(parse_float.(sline[6:8]))
            sum_spreads = parse_float(sline[9])
        end
        occursin(WOUT_MARKS.finalstate_end, line) && break
        line = readstrip(io)
    end
    return (; centers, spreads, sum_centers, sum_spreads)
end

"""
See [`_wout_parse_wf_center_spread`](@ref) for the format of the block to be parsed.
"""
@inline function _wout_parse_final_state(io::IO)
    return _wout_parse_wf_center_spread(io)
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
function _wout_parse_wannierize(io::IO)
    iter = Int[]
    centers = Vector{Vector{Vec3{Float64}}}()
    spreads = Vector{Vector{Float64}}()
    sum_centers = Vector{Vec3{Float64}}()
    sum_spreads = Float64[]
    ΩD = Float64[]
    ΩOD = Float64[]
    Ωtotal = Float64[]

    function append_iter!(idx::Int)
        c, s, sc, ss = _wout_parse_wf_center_spread(io)
        push!(centers, c)
        push!(spreads, s)
        push!(sum_centers, sc)
        push!(sum_spreads, ss)
        push!(iter, idx)
        sprd_line = ""
        while !eof(io)
            sprd_line = readstrip(io)
            occursin("<-- SPRD", sprd_line) && break
        end
        occursin("<-- SPRD", sprd_line) || error("No `<-- SPRD` found")
        sline = split(sprd_line)
        push!(ΩD, parse_float(sline[2]))
        push!(ΩOD, parse_float(sline[4]))
        push!(Ωtotal, parse_float(sline[6]))
        return nothing
    end

    while !eof(io)
        line = readstrip(io)
        if occursin(WOUT_MARKS.maxloc_end, line)
            # If not converged, w90 won't print
            # `<<< Wannierisation convergence criteria satisfied >>>`,
            # instead, it directly print "Final State".
            # Therefore, we use "Final State" as the marker
            # for the end of max loc iterations.
            break
        elseif occursin("Initial State", line)
            append_iter!(0)
        elseif occursin("Cycle:", line)
            append_iter!(parse(Int, split(line)[2]))
        end
    end
    return OrderedDict{String, Any}(
        "iter" => iter,
        "centers" => centers,
        "spreads" => spreads,
        "sum_centers" => sum_centers,
        "sum_spreads" => sum_spreads,
        "ΩD" => ΩD,
        "ΩOD" => ΩOD,
        "Ωtotal" => Ωtotal,
    )
end
