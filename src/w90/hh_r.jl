export write_HH_R

"""
Container for `prefix_HH_R.dat` data.

$(TYPEDEF)

# Fields

$(FIELDS)
"""
struct HHRDat{T<:Real,IT<:Integer}
    "Header line"
    header::String

    "`R` vectors of length `n_rvecs`"
    Rvectors::Vector{Vec3{IT}}

    "Degeneracy vector for each R vector, or `nothing`"
    Rdegens::Union{Vector{IT},Nothing}

    "Hamiltonian of length `n_rvecs`, each element is a matrix with shape `(n_wann, n_wann)`"
    H::Vector{Matrix{Complex{T}}}
end

function Base.show(io::IO, hhr::HHRDat)
    n_wann = isempty(hhr.H) ? 0 : size(hhr.H[1], 1)
    print(io, "HHRDat(n_Rvecs=$(length(hhr.H)), n_wann=$(n_wann))")
end

function Base.show(io::IO, ::MIME"text/plain", hhr::HHRDat)
    n_Rvecs = length(hhr.H)
    n_wann = isempty(hhr.H) ? 0 : size(hhr.H[1], 1)
    degen_str = if isnothing(hhr.Rdegens)
        "none"
    else
        degen_min = minimum(hhr.Rdegens)
        degen_max = maximum(hhr.Rdegens)
        "[$(degen_min), ..., $(degen_max)]"
    end

    print(
        io,
        """HHRDat(
          header: $(hhr.header)
          n_Rvecs: $(n_Rvecs)
          n_wann: $(n_wann)
          Rdegens: $(degen_str)
          H: Vector{Matrix{Complex}}($(n_wann)×$(n_wann))
        )""",
    )
end

"""
    $(SIGNATURES)

Write the real space Hamiltonian to a `prefix_HH_R.dat` file.

# Arguments
- `file`: The name of the output file, or an `IO`.
- `hhr`: a [`HHRDat`](@ref) struct

!!! note

    `Wannier90` `postw90.x` has a hidden input parameter `effective_model`,
    setting it to `true` and `postw90.x` will read this `HH_R.dat` to fill the
    real space Hamiltonian, and do subsequent Wannier interpolation, e.g.,
    in `BoltzWann`. However, the vanilla `postw90.x` code does not take into
    account the degeneracy of R vectors, and also does not use MDRS
    interpolation. I have modified the `postw90.x` code to use MDRS, and also
    changed a bit the number of digits for the Hamiltonian in `HH_R.dat`, so
    that it is the same as the `prefix_tb.dat` file, i.e., from Fortran
    `F12.6` to `E15.8`.

"""
function write_HH_R(io::IO, hhr::HHRDat)
    n_rvecs = length(hhr.Rvectors)
    n_rvecs > 0 || throw(ArgumentError("Rvectors must be non-empty"))
    length(hhr.H) == n_rvecs ||
        throw(DimensionMismatch("H and Rvectors must have the same length"))
    isnothing(hhr.Rdegens) ||
        length(hhr.Rdegens) == n_rvecs ||
        throw(DimensionMismatch("Rdegens and Rvectors must have the same length"))

    n_wann = size(hhr.H[1], 1)
    n_wann > 0 || throw(ArgumentError("H matrices must be non-empty"))
    all(size.(hhr.H) .== Ref((n_wann, n_wann))) ||
        throw(DimensionMismatch("H[ir] must all be n_wann * n_wann matrices"))

    vec2str(v) = join([@sprintf "%5d" x for x in v], "")

    write(io, strip(hhr.header), "\n")

    @printf(io, "%d\n", n_wann)
    @printf(io, "%d\n", n_rvecs)

    for ir in 1:n_rvecs
        for j in 1:n_wann
            for i in 1:n_wann
                h = hhr.H[ir][i, j]
                # 12.6f is the wannier90 default, however, I change it to
                # 15.8e so that it has the same accuracy as tb.dat file.
                @printf(
                    io,
                    "%s   %15.8e %15.8e\n",
                    vec2str([hhr.Rvectors[ir]..., i, j]),
                    real(h),
                    imag(h)
                )
            end
        end
    end

    return nothing
end

function write_HH_R(filename::AbstractString, hhr::HHRDat)
    n_rvecs = length(hhr.Rvectors)
    n_rvecs > 0 || throw(ArgumentError("Rvectors must be non-empty"))
    length(hhr.H) == n_rvecs ||
        throw(DimensionMismatch("H and Rvectors must have the same length"))
    isnothing(hhr.Rdegens) ||
        length(hhr.Rdegens) == n_rvecs ||
        throw(DimensionMismatch("Rdegens and Rvectors must have the same length"))

    n_wann = size(hhr.H[1], 1)
    n_wann > 0 || throw(ArgumentError("H matrices must be non-empty"))
    all(size.(hhr.H) .== Ref((n_wann, n_wann))) ||
        throw(DimensionMismatch("H[ir] must all be n_wann * n_wann matrices"))

    open(filename, "w") do io
        write_HH_R(io, hhr)
    end

    # the vanilla wannier90 code does not read the N array (the degeneracy
    # of R vector), and assume that N = 1 for all the R vectors.
    # I write it, and also implement the reading of wsvec.dat as well,
    # to use MDRS interpolation.
    if !isnothing(hhr.Rdegens)
        ndegen_filename = filename * ".ndegen"
        vec2str(v) = join([@sprintf "%5d" x for x in v], "")
        open(ndegen_filename, "w") do io
            # for aesthetic purpose, I write the N array in 15 columns
            n_col = 15  # 15 numbers per row
            for i in 0:(n_rvecs ÷ n_col - 1)
                s = i * n_col + 1  # start
                e = (i + 1) * n_col  # end
                @printf(io, "%s\n", vec2str(hhr.Rdegens[s:e]))
            end
            if (n_rvecs % n_col) > 0
                s = n_rvecs - n_rvecs % n_col + 1 # start
                e = n_rvecs  # end
                @printf(io, "%s\n", vec2str(hhr.Rdegens[s:e]))
            end
        end
    end
end
