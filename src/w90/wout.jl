using Printf: @printf

export read_wout

"""
    read_wout(filename::AbstractString)

Parse `wout` file.

Return lattice in Å (each column is a lattice vector),
atom positions in fractional coordinates (each column is a coordinate),
WF centers in Å, and WF spreads in Å^2.
"""
function read_wout(filename::AbstractString)
    io = open(filename)

    start_cell = "Lattice Vectors ("
    start_atom = "|   Site       Fractional Coordinate          Cartesian Coordinate"
    end_atom = "*----------------------------------------------------------------------------*"
    start_finalstate = "Final State"
    end_finalstate = "Sum of centres and spreads"

    ang_unit = false
    lattice = nothing
    atom_labels = nothing
    atom_positions = nothing
    centers = nothing
    spreads = nothing

    while !eof(io)
        line = strip(readline(io))

        if occursin("|  Length Unit", line)
            @assert occursin("Ang", line)
            ang_unit = true
            continue
        end

        if occursin(start_cell, line)
            @assert occursin("Ang", line)
            lattice = zeros(Float64, 3, 3)

            line = split(strip(readline(io)))
            @assert line[1] == "a_1"
            lattice[:, 1] = parse.(Float64, line[2:end])

            line = split(strip(readline(io)))
            @assert line[1] == "a_2"
            lattice[:, 2] = parse.(Float64, line[2:end])

            line = split(strip(readline(io)))
            @assert line[1] == "a_3"
            lattice[:, 3] = parse.(Float64, line[2:end])

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
            atom_positions = zeros(Float64, 3, n_atom)
            for (i, line) in enumerate(lines)
                line = split(line)
                @assert line[1] == "|" line
                push!(atom_labels, line[2])
                # cartesian
                # atom_positions[:, i] = parse.(Float64, line[8:10])
                # fractional
                atom_positions[:, i] = parse.(Float64, line[4:6])
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
            centers = zeros(Float64, 3, n_wann)
            spreads = zeros(Float64, n_wann)
            for (i, line) in enumerate(lines)
                line = split(line)
                @assert join(line[1:4], " ") == "WF centre and spread"
                @assert i == parse(Int, line[5])

                x = parse(Float64, replace(line[7], "," => ""))
                y = parse(Float64, replace(line[8], "," => ""))
                z = parse(Float64, replace(line[9], "," => ""))
                s = parse(Float64, line[11])

                centers[:, i] = [x, y, z]
                spreads[i] = s
            end

            continue
        end
    end

    close(io)

    @assert ang_unit

    return (; lattice, atom_labels, atom_positions, centers, spreads)
end
