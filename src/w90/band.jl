using Printf: @printf, @sprintf
using DelimitedFiles: readdlm

export read_w90_band, write_w90_band

function read_w90_band_kpt(filename::AbstractString)
    # in fractional coordinates
    kpoints = readdlm(filename, Float64; skipstart=1)
    weights = kpoints[:, 4]
    # remove weights, then transpose, last idx is kpt
    kpoints = Matrix(transpose(kpoints[:, 1:3]))
    return kpoints, weights
end

function read_w90_band_dat(filename::AbstractString)
    # 1st read to get n_kpts
    io = open(filename)
    n_kpts = 0
    while true
        line = strip(readline(io))
        isempty(line) && break
        n_kpts += 1
    end

    seekstart(io)
    dat = readdlm(io, Float64)
    close(io)

    x = reshape(dat[:, 1], n_kpts, :)[:, 1]
    E = reshape(dat[:, 2], n_kpts, :)
    # first dim is n_bands, 2nd dim is n_kpts
    E = Matrix(transpose(E))
    return x, E
end

function read_w90_band_labelinfo(filename::AbstractString)
    io = open(filename)
    labels = readlines(io)
    close(io)

    n_symm = length(labels)
    symm_idx = Vector{Int}(undef, n_symm)
    symm_label = Vector{String}(undef, n_symm)
    for (i, line) in enumerate(labels)
        lab, idx = split(line)[1:2]
        symm_idx[i] = parse(Int, idx)
        symm_label[i] = lab
    end

    return symm_idx, symm_label
end

"""
    read_w90_band(seedname::AbstractString)

Read `SEEDNAME_band.dat`, `SEEDNAME_band.kpt`, `SEEDNAME_band.labelinfo.dat`.

See also [`read_w90_band(seedname::AbstractString, recip_lattice::AbstractMatrix)`]
(@ref read_w90_band(seedname::AbstractString, recip_lattice::AbstractMatrix)).
"""
function read_w90_band(seedname::AbstractString)
    kpoints, _ = read_w90_band_kpt("$(seedname)_band.kpt")
    x, E = read_w90_band_dat("$(seedname)_band.dat")
    symm_idx, symm_label = read_w90_band_labelinfo("$(seedname)_band.labelinfo.dat")
    return (; kpoints, E, x, symm_idx, symm_label)
end

function write_w90_band_kpt(
    filename::AbstractString,
    kpoints::AbstractMatrix{T},
    weights::Union{Nothing,AbstractVector{T}}=nothing,
) where {T<:Real}
    n_kpts = size(kpoints, 2)
    size(kpoints, 1) == 3 || error("kpoints must be 3 x n_kpts")

    isnothing(weights) && (weights = ones(T, n_kpts))
    length(weights) == n_kpts || error("weights must be n_kpts")

    open(filename, "w") do io
        @printf(io, "       %5d\n", n_kpts)
        for ik in 1:n_kpts
            k = kpoints[:, ik]
            w = weights[ik]
            @printf(io, "  %10.6f  %10.6f  %10.6f   %10.6f\n", k..., w)
        end
    end
    @info "Written to $filename"
end

function write_w90_band_dat(
    filename::AbstractString, x::AbstractVector{T}, E::AbstractMatrix{T}
) where {T<:Real}
    n_bands, n_kpts = size(E)
    length(x) == n_kpts || error("x must be n_kpts")

    open(filename, "w") do io
        for ib in 1:n_bands
            for ik in 1:n_kpts
                @printf(io, " %15.8E %15.8E\n", x[ik], E[ib, ik])
            end
            @printf(io, "\n")
        end
    end
    @info "Written to $filename"
end

function write_w90_band_labelinfo(
    filename::AbstractString,
    symm_idx::AbstractVector{T},
    symm_label::AbstractVector{R},
    x::AbstractVector{S},
    kpoints::AbstractMatrix{S},
) where {T<:Integer,R<:AbstractString,S<:Real}
    n_symm = length(symm_idx)
    n_symm == length(symm_label) || error("symm_idx and symm_label must be same length")

    open(filename, "w") do io
        for i in 1:n_symm
            idx = symm_idx[i]
            @printf(
                io,
                "%2s %31d %20.10f %17.10f %17.10f %17.10f\n",
                symm_label[i],
                idx,
                x[idx],
                kpoints[:, idx]...
            )
        end
    end
    @info "Written to $filename"
end

"""
    write_w90_band(seedname, kpoints, E, x, symm_idx, symm_label)

Write `SEEDNAME_band.dat, SEEDNAME_band.kpt, SEEDNAME_band.labelinfo.dat`.

See also [`write_w90_band(seedname, kpi, E)`](@ref write_w90_band(seedname, kpi, E)).
"""
function write_w90_band(
    seedname::AbstractString,
    kpoints::AbstractMatrix{T},
    E::AbstractMatrix{T},
    x::AbstractVector{T},
    symm_idx::AbstractVector{R},
    symm_label::AbstractVector{S},
) where {T<:Real,R<:Integer,S<:AbstractString}
    size(kpoints, 2) == size(E, 2) || error("kpoints and E have different n_kpts")

    filename = "$(seedname)_band.kpt"
    write_w90_band_kpt(filename, kpoints)

    filename = "$(seedname)_band.dat"
    write_w90_band_dat(filename, x, E)

    filename = "$(seedname)_band.labelinfo.dat"
    write_w90_band_labelinfo(filename, symm_idx, symm_label, x, kpoints)

    println()
    return nothing
end
