# using GarishPrint: GarishPrint
# using Configurations: Configurations

# Base.@kwdef mutable struct AtomicWavefunction
#     # to differentiate same kind of atoms
#     atom_index::Int
#     # e.g. "Si", use "Si1", "Si2" to differentiate same kind but different types (spin up, down)
#     atom_label::String
#     # orbital label, e.g. "3S"
#     wfc_label::String

#     # quantum numbers
#     n::Int
#     l::Int
#     m::Int
# end

# Base.@kwdef mutable struct Projectabilities
#     n_kpts::Int
#     n_bands::Int
#     # number of atomic wavefunctions
#     num_wfcs::Int

#     # atomic wavefunction types, size: num_wfcs
#     wfcs_type::Vector{AtomicWavefunction}

#     # projectability data, size: n_kpts * n_bands * num_wfcs
#     proj::Array{Float64,3}
# end


# function read_qe_projwfcup(filename::String)
#     fdat = open(filename)

#     splitline() = split(strip(readline(fdat)))

#     # header
#     title = strip(readline(fdat), '\n')
#     nr1x, nr2x, nr3x, nr1, nr2, nr3, nat, ntyp = parse.(Int, splitline())
#     line = splitline()
#     ibrav = parse(Int, line[1])
#     celldm = parse.(Float64, line[2:end])
#     line = splitline()
#     # some version of projwfc.x output the unit_cell
#     if length(line) == 3
#         readline(fdat)
#         readline(fdat)
#         line = splitline()
#     end
#     gcutm, dual, ecutwfc = parse.(Float64, line[1:(end - 1)])
#     magicnum = parse(Int, line[end])
#     @assert magicnum == 9
#     atm = Vector{String}(undef, ntyp)
#     zv = zeros(Float64, ntyp)
#     for i in 1:ntyp
#         line = splitline()
#         nt = parse(Int, line[1])
#         @assert nt == i
#         atm[i] = line[2]
#         zv[i] = parse(Float64, line[3])
#     end
#     tau = zeros(Float64, 3, nat)
#     ityp = zeros(Int, nat)
#     for i in 1:nat
#         line = splitline()
#         na = parse(Int, line[1])
#         @assert na == i
#         tau[:, i] = parse.(Float64, line[2:4])
#         ityp[i] = parse(Int, line[5])
#     end
#     natomwfc, nkstot, nbnd = parse.(Int, splitline())
#     parsebool(s::Union{String,SubString}) = lowercase(s) == "t" ? true : false
#     noncolin, lspinorb = parsebool.(splitline())
#     @assert !noncolin && !lspinorb

#     # projection data
#     nlmchi = Vector{Dict}()
#     proj = zeros(Float64, nkstot, nbnd, natomwfc)
#     for iw in 1:natomwfc
#         line = splitline()
#         nwfc = parse(Int, line[1])
#         @assert nwfc == iw
#         na = parse(Int, line[2])
#         atm_name = line[3]
#         @assert atm_name == atm[ityp[na]]
#         els = line[4]
#         n, l, m = parse.(Int, line[5:end])
#         push!(nlmchi, Dict("na" => na, "els" => els, "n" => n, "l" => l, "m" => m))
#         for ik in 1:nkstot
#             for ib in 1:nbnd
#                 line = splitline()
#                 k, b = parse.(Int, line[1:2])
#                 @assert k == ik && b == ib
#                 p = parse(Float64, line[3])
#                 proj[ik, ib, iw] = p
#             end
#         end
#     end

#     wfcs_type = Vector{AtomicWavefunction}(undef, natomwfc)
#     for i in 1:natomwfc
#         atom_index = nlmchi[i]["na"]
#         atom_label = atm[ityp[atom_index]]
#         wfc_label = nlmchi[i]["els"]
#         n = nlmchi[i]["n"]
#         l = nlmchi[i]["l"]
#         m = nlmchi[i]["m"]
#         wfc = AtomicWavefunction(atom_index, atom_label, wfc_label, n, l, m)
#         wfcs_type[i] = wfc
#     end

#     return Projectabilities(nkstot, nbnd, natomwfc, wfcs_type, proj)
# end
