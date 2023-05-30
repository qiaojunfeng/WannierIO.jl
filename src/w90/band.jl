using Printf: @printf, @sprintf
using DelimitedFiles: readdlm

export read_w90_band, write_w90_band

"""
    read_w90_band_kpt(filename)

Read `seedname_band.kpt` file.

# Return
- `kpoints`: a vector of length `n_kpts`, fractional coordinates
- `weights`: a vector of length `n_kpts`, weights of kpoints
"""
function read_w90_band_kpt(filename::AbstractString)
    # in fractional coordinates
    kpoints = readdlm(filename, Float64; skipstart=1)
    weights = kpoints[:, 4]
    # remove weights
    kpoints = map(1:size(kpoints, 1)) do i
        Vec3(kpoints[i, 1:3])
    end
    return kpoints, weights
end

"""
    read_w90_band_dat(filename)

Read `seedname_band.dat` file.

# Return
- `x`: `n_kpts`, x axis value, in cartesian length
- `E`: length-`n_kpts` vector, each elemnt is a length-`n_bands` vector of band energies
"""
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
    E = [Vector(e) for e in eachrow(E)]
    return x, E
end

"""
    read_w90_band_labelinfo(filename::AbstractString)

Read `seedname_band.labelinfo` file.

# Return
- `symm_idx`: index of high-symmetry points in `seedname_band.dat`
- `symm_label`: name of high-symmetry points
"""
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

# Return
- `kpoints`: `3 * n_kpts`, fractional coordinates
- `E`: length-`n_kpts` vector, each element is a length-`n_bands` vector of band energies
- `x`: `n_kpts`, x axis value, in cartesian length
- `symm_idx`: index of high-symmetry points in `seedname_band.dat`
- `symm_label`: name of high-symmetry points
"""
function read_w90_band(seedname::AbstractString)
    kpoints, _ = read_w90_band_kpt("$(seedname)_band.kpt")
    x, E = read_w90_band_dat("$(seedname)_band.dat")
    symm_idx, symm_label = read_w90_band_labelinfo("$(seedname)_band.labelinfo.dat")
    return (; kpoints, E, x, symm_idx, symm_label)
end

"""
    write_w90_band_kpt(filename, kpoints, weights=nothing)

Write `seedname_band.kpt` file.

# Arguments
- `filename`: filename of `seedname_band.kpt`
- `kpoints`: length-`n_kpts` vector, fractional coordinates
- `weights`: `n_kpts`, optional, weights of kpoints
"""
function write_w90_band_kpt(
    filename::AbstractString,
    kpoints::AbstractVector{<:Vec3{<:Real}},
    weights::Union{Nothing,AbstractVector{<:Real}}=nothing,
)
    n_kpts = length(kpoints)

    isnothing(weights) && (weights = ones(n_kpts))
    length(weights) == n_kpts || error("weights must be n_kpts")

    open(filename, "w") do io
        @printf(io, "       %5d\n", n_kpts)
        for ik in 1:n_kpts
            k = kpoints[ik]
            w = weights[ik]
            @printf(io, "  %10.6f  %10.6f  %10.6f   %10.6f\n", k..., w)
        end
    end
    @info "Written to $filename"
end

"""
    write_w90_band_dat(filename, x, E)

Write `seedname_band.dat` file.

# Arguments
- `filename`: filename of `seedname_band.dat`
- `x`: `n_kpts`, x axis value, in cartesian length
- `E`: length-`n_kpts` vector, each element is a length-`n_bands` vector of band energies
"""
function write_w90_band_dat(
    filename::AbstractString,
    x::AbstractVector{<:Real},
    E::AbstractVector{<:AbstractVector{<:Real}},
)
    n_bands = length(E[1])
    n_kpts = length(E)
    length(x) == n_kpts || error("x must be n_kpts")

    open(filename, "w") do io
        for ib in 1:n_bands
            for ik in 1:n_kpts
                @printf(io, " %15.8E %15.8E\n", x[ik], E[ik][ib])
            end
            @printf(io, "\n")
        end
    end
    @info "Written to $filename"
end

"""
    write_w90_band_labelinfo(filename, symm_idx, symm_label, x, kpoints)

Write `seedname_band.labelinfo` file.

# Arguments
- `filename`: filename of `seedname_band.labelinfo`
- `symm_idx`: index of high-symmetry points in `seedname_band.dat`
- `symm_label`: name of high-symmetry points
- `x`: `n_kpts`, x axis value, in cartesian length
- `kpoints`: length-`n_kpts` vector, fractional coordinates
"""
function write_w90_band_labelinfo(
    filename::AbstractString,
    symm_idx::AbstractVector{<:Integer},
    symm_label::AbstractVector{<:AbstractString},
    x::AbstractVector{<:Real},
    kpoints::AbstractVector{<:Vec3{<:Real}},
)
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
                kpoints[idx]...
            )
        end
    end
    @info "Written to $filename"
end

"""
    write_w90_band(seedname, kpoints, E, x, symm_idx, symm_label)

Write `SEEDNAME_band.dat, SEEDNAME_band.kpt, SEEDNAME_band.labelinfo.dat`.

# Arguments
- `seedname`: seedname of `SEEDNAME_band.dat, SEEDNAME_band.kpt, SEEDNAME_band.labelinfo.dat`
- `kpoints`: length-`n_kpts` vector, fractional coordinates
- `E`: length-`n_kpts` vector, each element is a length-`n_bands` vector of band energies
- `x`: `n_kpts`, x axis value, in cartesian length
- `symm_idx`: index of high-symmetry points in `seedname_band.dat`
- `symm_label`: name of high-symmetry points
"""
function write_w90_band(
    seedname::AbstractString,
    kpoints::AbstractVector,
    E::AbstractVector,
    x::AbstractVector,
    symm_idx::AbstractVector,
    symm_label::AbstractVector,
)
    length(kpoints) == length(E) || error("kpoints and E have different n_kpts")

    filename = "$(seedname)_band.kpt"
    write_w90_band_kpt(filename, kpoints)

    filename = "$(seedname)_band.dat"
    write_w90_band_dat(filename, x, E)

    filename = "$(seedname)_band.labelinfo.dat"
    write_w90_band_labelinfo(filename, symm_idx, symm_label, x, kpoints)

    println()
    return nothing
end
