using Printf: @printf, @sprintf
using DelimitedFiles: readdlm

export read_w90_band, write_w90_band

"""
    read_w90_band_kpt(filename)

Read a `prefix_band.kpt` file.

# Return
- `kpoints`: a vector of length `n_kpts`, fractional coordinates
- `kweights`: a vector of length `n_kpts`, weights of kpoints
"""
function read_w90_band_kpt(filename::AbstractString)
    # in fractional coordinates
    kpoints = readdlm(filename, Float64; skipstart=1)
    kweights = kpoints[:, 4]
    # remove weights
    kpoints = map(1:size(kpoints, 1)) do i
        Vec3(kpoints[i, 1:3])
    end
    return (; kpoints, kweights)
end

"""
    read_w90_band_dat(filename)

Read `prefix_band.dat` file.

# Return
- `x`: `n_kpts`, x axis value of kpath, in cartesian length
- `eigenvalues`: length-`n_kpts` vector, each elemnt is a length-`n_bands` vector of band energies
"""
function read_w90_band_dat(filename::AbstractString)
    res = open(filename) do io
        # Unfortunately I need to read the whole file twice:
        # 1st time to get n_kpts, 2nd time to get data
        n_kpts = 0
        while true
            line = strip(readline(io))
            isempty(line) && break
            n_kpts += 1
        end

        seekstart(io)
        dat = readdlm(io, Float64)

        x = reshape(dat[:, 1], n_kpts, :)[:, 1]
        eigenvalues = reshape(dat[:, 2], n_kpts, :)
        eigenvalues = [Vector(e) for e in eachrow(eigenvalues)]
        return (; x, eigenvalues)
    end
    return res
end

"""
    read_w90_band_labelinfo(filename)

Read `prefix_band.labelinfo` file.

# Return
- `symm_point_indices`: index of high-symmetry points in `prefix_band.dat`
- `symm_point_labels`: name of high-symmetry points
"""
function read_w90_band_labelinfo(filename::AbstractString)
    labels = open(filename) do io
        readlines(io)
    end

    n_symm = length(labels)
    symm_point_indices = Vector{Int}(undef, n_symm)
    symm_point_labels = Vector{String}(undef, n_symm)
    for (i, line) in enumerate(labels)
        lab, idx = split(line)[1:2]
        symm_point_indices[i] = parse(Int, idx)
        symm_point_labels[i] = lab
    end

    return (; symm_point_indices, symm_point_labels)
end

"""
    read_w90_band(prefix)

Read `prefix_band.dat`, `prefix_band.kpt`, `prefix_band.labelinfo.dat`.

# Arguments
- `prefix`: *prefix* of the filenames (or called seedname in wannier90), NOT the full filename.

# Return
- `x`: `n_kpts`, x axis value of kpath, in cartesian length
- `eigenvalues`: length-`n_kpts` vector, each element is a length-`n_bands` vector of band energies
- `kpoints`: a vector of length `n_kpts`, fractional coordinates
- `kweights`: a vector of length `n_kpts`, weights of kpoints
- `symm_point_indices`: index of high-symmetry points in `prefix_band.dat`
- `symm_point_labels`: name of high-symmetry points
"""
function read_w90_band(prefix::AbstractString)
    band_dat = prefix * "_band.dat"
    band_kpt = prefix * "_band.kpt"
    band_labelinfo = prefix * "_band.labelinfo.dat"

    dat = read_w90_band_dat(band_dat)
    kpt = read_w90_band_kpt(band_kpt)
    labelinfo = read_w90_band_labelinfo(band_labelinfo)

    n_kpts = length(kpt.kpoints)
    n_symm = length(labelinfo.symm_point_indices)
    @info "Reading Wannier90 band files" band_dat band_kpt band_labelinfo n_kpts n_symm

    return (; dat..., kpt..., labelinfo...)
end

"""
Wannier90 default kweights in `prefix_band.kpt` is all 1.0.
"""
default_band_kpt_kweights(kpoints::AbstractVector) = ones(length(kpoints))

"""
    write_w90_band_kpt(filename; kpoints, kweights=default_band_kpt_kweights(kpoints))

Write `prefix_band.kpt` file.

# Arguments
- `filename`: filename of `prefix_band.kpt`

# Keyword Arguments
- `kpoints`: length-`n_kpts` vector, fractional coordinates
- `kweights`: `n_kpts`, optional, weights of kpoints
"""
function write_w90_band_kpt(
    filename::AbstractString;
    kpoints::AbstractVector,
    kweights::AbstractVector=default_band_kpt_kweights(kpoints),
)
    n_kpts = length(kpoints)
    length(kweights) == n_kpts || error("kweights must have same length as kpoints")

    open(filename, "w") do io
        @printf(io, "       %5d\n", n_kpts)
        for (k, w) in zip(kpoints, kweights)
            length(k) == 3 || error("kpoint must be 3-vector")
            @printf(io, "  %10.6f  %10.6f  %10.6f   %10.6f\n", k..., w)
        end
    end
end

"""
    write_w90_band_dat(filename; x, eigenvalues)

Write `prefix_band.dat` file.

# Arguments
- `filename`: filename of `prefix_band.dat`

# Keyword Arguments
- `x`: `n_kpts`, x axis value, in cartesian length
- `eigenvalues`: length-`n_kpts` vector, each element is a length-`n_bands` vector of band energies
"""
function write_w90_band_dat(
    filename::AbstractString;
    x::AbstractVector,
    eigenvalues::AbstractVector{<:AbstractVector},
)
    n_kpts = length(eigenvalues)
    @assert n_kpts > 0 "eigenvalues is empty"
    n_bands = length(eigenvalues[1])
    length(x) == n_kpts || error("x must has same length as eigenvalues")

    open(filename, "w") do io
        for ib in 1:n_bands
            for ik in 1:n_kpts
                @printf(io, " %15.8E %15.8E\n", x[ik], eigenvalues[ik][ib])
            end
            @printf(io, "\n")
        end
    end
end

"""
    write_w90_band_labelinfo(filename; x, kpoints, symm_point_indices, symm_point_labels)

Write `prefix_band.labelinfo` file.

# Arguments
- `filename`: filename of `prefix_band.labelinfo`

# Keyword Arguments
- `x`: `n_kpts`-vector, x axis value, in cartesian length
- `kpoints`: length-`n_kpts` vector, fractional coordinates
- `symm_point_indices`: index of high-symmetry points in `prefix_band.dat`
- `symm_point_labels`: name of high-symmetry points
"""
function write_w90_band_labelinfo(
    filename::AbstractString;
    x::AbstractVector{<:Real},
    kpoints::AbstractVector,
    symm_point_indices::AbstractVector{<:Integer},
    symm_point_labels::AbstractVector,
)
    n_symm = length(symm_point_indices)
    n_symm == length(symm_point_labels) ||
        error("symm_idx and symm_label must have same length")

    open(filename, "w") do io
        for i in 1:n_symm
            idx = symm_point_indices[i]
            kpt = kpoints[idx]
            length(kpt) == 3 || error("kpoint must be 3-vector")
            @printf(
                io,
                "%2s %31d %20.10f %17.10f %17.10f %17.10f\n",
                symm_point_labels[i],
                idx,
                x[idx],
                kpt...
            )
        end
    end
end

"""
    write_w90_band(prefix; x, eigenvalues, kpoints, kweights, symm_point_indices, symm_point_labels)

Write `prefix_band.dat, prefix_band.kpt, prefix_band.labelinfo.dat`.

# Arguments
- `prefix`: prefix of `prefix_band.dat, prefix_band.kpt, prefix_band.labelinfo.dat`

# Keyword Arguments
- `x`: `n_kpts`, x axis value, in cartesian length
- `eigenvalues`: length-`n_kpts` vector, each element is a length-`n_bands` vector of band energies
- `kpoints`: length-`n_kpts` vector, fractional coordinates
- `kweights`: a vector of length `n_kpts`, weights of kpoints
- `symm_point_indices`: index of high-symmetry points in `prefix_band.dat`
- `symm_point_labels`: name of high-symmetry points
"""
function write_w90_band(
    prefix::AbstractString;
    x::AbstractVector,
    eigenvalues::AbstractVector{<:AbstractVector},
    kpoints::AbstractVector,
    kweights::AbstractVector=default_band_kpt_kweights(kpoints),
    symm_point_indices::AbstractVector,
    symm_point_labels::AbstractVector,
)
    n_kpts = length(kpoints)
    length(eigenvalues) == n_kpts || error("kpoints and eigenvalues have different n_kpts")
    length(kweights) == n_kpts || error("kpoints and kweights have different n_kpts")
    length(x) == n_kpts || error("kpoints and x have different n_kpts")
    n_symm = length(symm_point_indices)
    n_symm == length(symm_point_labels) ||
        error("symm_idx and symm_label must have same length")

    band_kpt = "$(prefix)_band.kpt"
    band_dat = "$(prefix)_band.dat"
    band_labelinfo = "$(prefix)_band.labelinfo.dat"

    @info "Writing Wannier90 band files" band_kpt band_dat band_labelinfo n_kpts n_symm

    write_w90_band_dat(band_dat; x, eigenvalues)
    write_w90_band_kpt(band_kpt; kpoints, kweights)
    return write_w90_band_labelinfo(
        band_labelinfo; x, kpoints, symm_point_indices, symm_point_labels
    )
end
