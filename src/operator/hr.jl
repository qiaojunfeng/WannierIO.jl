export read_w90_hr, write_w90_hr

"""
    read_w90_hr(file)

Read Wannier90 `prefix_hr.dat` (and optional matching `prefix_wsvec.dat`).

The Rvectors and degeneracies will be reduced by [`RvectorReducer`](@ref)
automatically, resulting in a [`OperatorPack`](@ref) with just the Rvectors
and the `H` operator matrices defined on them.

If a corresponding `wsvec.dat` file is found, the MDRS Tvectors and degeneracies
will also be reduced automatically.

Because `prefix_hr.dat` does not store lattice vectors, the returned
[`OperatorPack`](@ref) uses a `3×3` lattice filled with `NaN` as a sentinel for
"lattice unavailable".

As a consequence, the Wannier interpolation will be just a simple inverse
Fourier transform without any complications from R-degeneracies or MDRS
T vectors and degeneracies.

See also [`read_w90_hr_dat`](@ref) and [`read_w90_wsvec`](@ref) for reading just
the `prefix_hr.dat` and `prefix_wsvec.dat` files without any reduction.

This high-level API is intentionally for native Wannier90 text files only.
For HDF5/JLD2/Zarr backends, use [`read_operator`](@ref).
"""
function read_w90_hr(file::AbstractString)
    hrdat = read_w90_hr_dat(file)
    wsvec_path = _wsvec_path_from_hr(file)
    if isfile(wsvec_path)
        wsvec = read_w90_wsvec(wsvec_path)
        return pack(hrdat, wsvec)
    else
        @warn("No corresponding wsvec file found for hr file at $file.")
        return pack(hrdat)
    end
end

"""
    write_w90_hr(file, hrdat)
    write_w90_hr(file, hrdat, wsvec)

Write Wannier90 Hamiltonian text data.

This high-level API is intentionally for native Wannier90 text files only.
For HDF5/JLD2/Zarr backends, use [`write_operator`](@ref).

For [`W90Dat`](@ref):
- `write_w90_hr(file, hrdat)` writes only `prefix_hr.dat`.
- `write_w90_hr(file, hrdat, wsvec)` writes both `prefix_hr.dat` and `prefix_wsvec.dat`.
"""
function write_w90_hr end

function write_w90_hr(file::AbstractString, hrdat::HrDat)
    write_w90_hr_dat(file, hrdat)
    return nothing
end

function write_w90_hr(file::AbstractString, hrdat::HrDat, wsvec::WsvecDat)
    write_w90_hr_dat(file, hrdat)
    write_w90_wsvec(_wsvec_path_from_hr(file), wsvec)
    return nothing
end

function _wsvec_path_from_hr(hrpath::AbstractString)
    if endswith(hrpath, "_hr.dat")
        return hrpath[1:(end - 7)] * "_wsvec.dat"
    else
        error("Not a valid hr.dat filename: $hrpath")
    end
end

function write_w90_hr(file::AbstractString, pack::OperatorPack)
    hrdat = HrDat(pack)
    wsvec = WsvecDat(pack)
    write_w90_hr(file, hrdat, wsvec)
    return nothing
end

function _missing_lattice(::Type{T}) where {T <: Real}
    return Mat3{T}(fill(T(NaN), 3, 3))
end

function RvectorReducer(hrdat::HrDat)
    return RvectorReducer(hrdat.Rvectors, hrdat.Rdegens)
end

function RvectorReducer(hrdat::HrDat, wsvec::WsvecDat)
    hrdat.Rvectors == wsvec.Rvectors ||
        error("Rvectors in hr.dat and wsvec.dat are inconsistent")

    if wsvec.mdrs
        return RvectorReducer(hrdat.Rvectors, hrdat.Rdegens, wsvec.Tvectors, wsvec.Tdegens)
    end
    return RvectorReducer(hrdat.Rvectors, hrdat.Rdegens)
end

function pack(hrdat::HrDat)
    reducer = RvectorReducer(hrdat)
    return pack(hrdat, reducer)
end

function pack(hrdat::HrDat, wsvec::WsvecDat)
    reducer = RvectorReducer(hrdat, wsvec)
    return pack(hrdat, reducer)
end

function pack(hrdat::HrDat{T}, reducer::AbstractRvectorReducer) where {T <: Real}
    Tr = float(T)
    lattice = _missing_lattice(Tr)
    Rvectors = Vec3{Int}.(reducer.Rvectors)
    ops = OrderedDict("H" => [Matrix{Complex{Tr}}(M) for M in reducer(hrdat.H)])
    return OperatorPack(hrdat.header, lattice, Rvectors, ops)
end

function HrDat(pack::OperatorPack)
    haskey(pack.operators, "H") || error("missing operator `H`")
    H = pack.operators["H"]
    Tc = Complex{real(eltype(eltype(H)))}
    return HrDat(
        pack.header,
        Vec3{Int}.(pack.Rvectors),
        ones(Int, pack.n_Rvecs),
        [Matrix{Tc}(O) for O in H],
    )
end
