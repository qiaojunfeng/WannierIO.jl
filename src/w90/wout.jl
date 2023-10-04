export read_wout

"""
    $(SIGNATURES)

Parse wannire90 `wout` file.

# Return
- `lattice`: each column is a lattice vector in Å
- `recip_lattice`: each column is a reciprocal lattice vector in Å⁻¹
- `atom_labels`: atomic symbols
- `atom_positions`: in fractional coordinates
- `centers`: final each WF centers in Å
- `spreads`: final each WF spreads in Å²
- `ΩI`, `ΩD`, `ΩOD`, `Ωtotal`: final spread (components) in Å²
"""
function read_wout(filename::AbstractString)
    return open(filename) do io
        start_lattice = "Lattice Vectors ("
        start_recip = "Reciprocal-Space Vectors ("
        start_atom = "|   Site       Fractional Coordinate          Cartesian Coordinate"
        end_atom = "*----------------------------------------------------------------------------*"
        start_finalstate = "Final State"
        end_finalstate = "Sum of centres and spreads"
        marker_wfc = "WF centre and spread"
        marker_ΩI = "Omega I      ="
        marker_ΩD = "Omega D      ="
        marker_ΩOD = "Omega OD     ="
        marker_Ωtotal = "Omega Total  ="

        ang_unit = false
        lattice = nothing
        recip_lattice = nothing
        atom_labels = nothing
        atom_positions = nothing
        centers = nothing
        spreads = nothing
        ΩI = nothing
        ΩD = nothing
        ΩOD = nothing
        Ωtotal = nothing

        while !eof(io)
            line = strip(readline(io))

            if occursin(start_lattice, line)
                @assert occursin("Ang", line)
                ang_unit = true
                lattice = zeros(Float64, 3, 3)

                line = split(strip(readline(io)))
                @assert line[1] == "a_1"
                lattice[:, 1] = parse_float.(line[2:end])

                line = split(strip(readline(io)))
                @assert line[1] == "a_2"
                lattice[:, 2] = parse_float.(line[2:end])

                line = split(strip(readline(io)))
                @assert line[1] == "a_3"
                lattice[:, 3] = parse_float.(line[2:end])

                continue
            end

            if occursin(start_recip, line)
                @assert occursin("Ang^-1", line)
                recip_lattice = zeros(Float64, 3, 3)

                line = split(strip(readline(io)))
                @assert line[1] == "b_1"
                recip_lattice[:, 1] = parse_float.(line[2:end])

                line = split(strip(readline(io)))
                @assert line[1] == "b_2"
                recip_lattice[:, 2] = parse_float.(line[2:end])

                line = split(strip(readline(io)))
                @assert line[1] == "b_3"
                recip_lattice[:, 3] = parse_float.(line[2:end])

                continue
            end

            if occursin(start_atom, line)
                @assert occursin("Ang", line)
                readline(io)

                lines = Vector{String}()
                line = strip(readline(io))
                while line != end_atom
                    push!(lines, line)
                    line = strip(readline(io))
                end

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

                continue
            end

            if occursin(start_finalstate, line)
                lines = Vector{String}()
                line = strip(readline(io))
                while !occursin(end_finalstate, line)
                    push!(lines, line)
                    line = strip(readline(io))
                end

                n_wann = length(lines)
                centers = zeros(Vec3{Float64}, n_wann)
                spreads = zeros(Float64, n_wann)
                for (i, line) in enumerate(lines)
                    @assert startswith(line, marker_wfc) line
                    line = split(line, r"[,()]")

                    idx = strip(chopprefix(line[1], marker_wfc))
                    idx = parse(Int, idx)
                    @assert i == idx "Wrong WF index at $(line[1])"

                    vals = map(parse_float, line[2:end])
                    centers[i] = Vec3(vals[1:3])
                    spreads[i] = vals[4]
                end

                continue
            end

            if occursin(marker_ΩI, line)
                ΩI = parse_float(split(line, marker_ΩI)[2])
                continue
            end
            if occursin(marker_ΩD, line)
                ΩD = parse_float(split(line, marker_ΩD)[2])
                continue
            end
            if occursin(marker_ΩOD, line)
                ΩOD = parse_float(split(line, marker_ΩOD)[2])
                continue
            end
            if occursin(marker_Ωtotal, line)
                Ωtotal = parse_float(split(line, marker_Ωtotal)[2])
                continue
            end
        end

        @assert ang_unit "wout unit is not Angstrom, not supported yet"
        lattice = Mat3(lattice)
        recip_lattice = Mat3(recip_lattice)
        return (;
            lattice,
            recip_lattice,
            atom_labels,
            atom_positions,
            centers,
            spreads,
            ΩI,
            ΩD,
            ΩOD,
            Ωtotal,
        )
    end
end
