using LinearAlgebra

"""
    $(SIGNATURES)

Read Quantum ESPRESSO `bands.x` output data file.

The data file has format
```
 &plot nbnd=  20, nks=   380 /
           -0.500000  0.500000  0.500000
   -3.320   -0.666    5.173    5.173    7.994    9.725    9.725   14.147   16.993   16.993
   17.841   17.841   17.902   19.666   25.961   26.563   28.186   28.186   28.368   28.368
           -0.495000  0.495000  0.495000
   -3.322   -0.664    5.173    5.173    7.994    9.725    9.725   14.148   16.980   16.980
...
```
"""
function read_qe_band(filename::AbstractString)
    res = open(filename) do io
        line = readline(io)
        regex = r"&plot nbnd=\s*(\d+), nks=\s*(\d+) /"
        m = match(regex, line)
        if m !== nothing
            n_bands, n_kpts = parse.(Int, m.captures)
        else
            # this is my customized version, with `alat` added to header,
            # so we can return kpoints in Å⁻¹ unit instead of arbitrary
            regex = r"&plot nbnd=\s*(\d+), nks=\s*(\d+) alat=\s*([+-]?([0-9]*[.])?[0-9]+) /"
            m = match(regex, line)
            n_bands, n_kpts = parse.(Int, m.captures[1:2])
            alat = parse.(Float64, m.captures[3])
        end

        kpoints = Vec3{Float64}[]
        eigenvalues = Vector{Float64}[]

        for _ in 1:n_kpts
            # QE kpt are in cartesian coordinates, but scaled by `alat`
            kpt = parse.(Float64, split(readline(io)))
            push!(kpoints, kpt)
            ib = 1
            eig = zeros(Float64, n_bands)
            while ib <= n_bands
                e = parse.(Float64, split(readline(io)))
                n_e = length(e)
                eig[ib:(ib + n_e - 1)] = e
                ib += n_e
            end
            @assert ib == (n_bands + 1)
            push!(eigenvalues, eig)
        end
        @assert eof(io)

        return (; kpoints, eigenvalues)
    end

    return res
end

"""
    $(SIGNATURES)

Guess high symmetry points from kpoint coordinates.

If there is angle between two consecutive kpoints, then
it is labeled as a high-symmetry point.

# Arguments
- `kpoints`: Vector of `Vector` or `Vec3`, in Cartesian coordinates

# Keyword Arguments
- `atol`: Absolute tolerance for checking cross product of two vectors

# Returns
- `symm_point_indices`: Vector of indices of high-symmetry points
- `symm_point_labels`: Vector of labels of high-symmetry points, for the moment
  it is empty
"""
function guess_kpath(kpoints::AbstractVector{<:AbstractVector}; atol=2e-6)
    # of course index starts from 1
    symm_point_indices = Vector{Int}()
    symm_point_labels = Vector{String}()

    n_kpts = length(kpoints)
    if n_kpts == 0
        return (; symm_point_indices, symm_point_labels)
    end

    # push the first kpt
    push!(symm_point_indices, 1)
    push!(symm_point_labels, "")

    for ik in 2:(n_kpts - 1)
        u = kpoints[ik] - kpoints[ik - 1]
        v = kpoints[ik + 1] - kpoints[ik]
        if !all(isapprox.(cross(u, v), 0; atol))
            push!(symm_point_indices, ik)
            push!(symm_point_labels, "")
        end
    end

    # push the last kpt
    push!(symm_point_indices, n_kpts)
    push!(symm_point_labels, "")

    return (; symm_point_indices, symm_point_labels)
end
