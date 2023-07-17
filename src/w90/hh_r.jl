# TODO update this file, use consistent variable names

"""
    $(SIGNATURES)

Write the real space Hamiltonian to a `prefix_HH_R.dat` file.

# Arguments
- `filename`: usually `prefix_HH_R.dat`
- `H`: a `n_wann * n_wann * n_rvecs` array of Hamiltonian
- `R`: a `n_rvecs * 3` array of integers

# Keyword arguments
- `N`: a `n_rvecs` vector of integers, the degeneracy of each R vector
- `header`: a string, the header of the file

!!! note

    `Wanier90` `postw90.x` has a hidden input parameter `effective_model`,
    setting it to `true` and `postw90.x` will read this `HH_R.dat` to fill the
    real space Hamiltonian, and do subsequent Wannier interpolation, e.g.,
    in `BoltzWann`. However, the vanila `postw90.x` code does not take into
    account the degeneracy of R vectors, and also does not use MDRS
    interpolation. I have modified the `postw90.x` code to use MDRS, and also
    changed a bit the number of digits for the Hamiltonian in `HH_R.dat`, so
    that it is the same as the `prefix_tb.dat` file, i.e., from Fortran
    `F12.6` to `E15.8`.

"""
function write_HH_R(
    filename::AbstractString,
    H::AbstractArray{T,3},
    R::AbstractMatrix{IT};
    N::Union{AbstractVector{IT},Nothing}=nothing,
    header=default_header(),
) where {T<:Complex,IT<:Integer}
    n_wann, _, n_rvecs = size(H)
    size(H, 2) == n_wann || error("H must be a n_wann * n_wann * n_rvecs matrix")
    size(R) == (3, n_rvecs) || error("R must be a 3 * n_rvecs matrix")
    N === nothing || length(N) == n_rvecs || error("N must be a n_rvecs vector")

    io = open(filename, "w")

    write(io, header, "\n")

    @printf(io, "%d\n", n_wann)
    @printf(io, "%d\n", n_rvecs)

    vec2str(v) = join([@sprintf "%5d" x for x in v], "")

    for ir in 1:n_rvecs
        for j in 1:n_wann
            for i in 1:n_wann
                h = H[i, j, ir]
                # 12.6f is the wannier90 default, however, I change it to
                # 15.8e so that it has the same accuracy as tb.dat file.
                @printf(
                    io,
                    "%s   %15.8e %15.8e\n",
                    vec2str([R[:, ir]..., i, j]),
                    real(h),
                    imag(h)
                )
            end
        end
    end
    close(io)
    @info "Written to file: $(filename)"

    # the vanila wannier90 code does not read the N array (the degeneracy
    # of R vector), and assume that N = 1 for all the R vectors.
    # I write it, and also implement the reading of wsvec.dat as well,
    # to use MDRS interpolation.
    if N !== nothing
        filename *= ".ndegen"
        io = open(filename, "w")

        # for aesthetic purpose, I write the N array in 15 columns
        n_col = 15  # 15 numbers per row
        for i in 0:(n_rvecs ÷ n_col - 1)
            s = i * n_col + 1  # start
            e = (i + 1) * n_col  # end
            @printf(io, "%s\n", vec2str(N[s:e]))
        end
        if (n_rvecs % n_col) > 0
            s = n_rvecs - n_rvecs % n_col + 1 # start
            e = n_rvecs  # end
            @printf(io, "%s\n", vec2str(N[s:e]))
        end

        @info "Written to file: $(filename)"
        close(io)
    end
end
