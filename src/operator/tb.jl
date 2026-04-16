export read_w90_tb, write_w90_tb

"""
    read_w90_tb(file)

Read Wannier90 `prefix_tb.dat` (and optional matching `prefix_wsvec.dat`).

The Rvectors and degeneracies will be reduced by [`RvectorReducer`](@ref)
automatically, resulting in a [`OperatorPack`](@ref) with just the Rvectors
and the operator matrices defined on them.

If a corresponding `wsvec.dat` file is found, the MDRS Tvectors and degeneracies
will also be reduced automatically.

As a consequence, the Wannier interpolation will be just a simple inverse
Fourier transform without any complications from R-degeneracies or MDRS
T vectors and degeneracies.

See also [`read_w90_tb_dat`](@ref) and [`read_w90_wsvec_dat`](@ref) for reading just
the `prefix_tb.dat` and `prefix_wsvec.dat` files without any reduction.

This high-level API is intentionally for native Wannier90 text files only.
For HDF5/JLD2/Zarr backends, use [`read_operator`](@ref).
"""
function read_w90_tb end

function read_w90_tb(file::AbstractString)
    tbdat = read_w90_tb_dat(file)
    wsvec_path = _wsvec_path_from_tb(file)
    if isfile(wsvec_path)
        wsvec = read_w90_wsvec_dat(wsvec_path)
        return pack(tbdat, wsvec)
    else
        @warn("No corresponding wsvec file found for tb file at $file.")
        return pack(tbdat)
    end
end

"""
    write_w90_tb(file, tbdat)
    write_w90_tb(file, tbdat, wsvec)

Write Wannier90 tight-binding text data.

This high-level API is intentionally for native Wannier90 text files only.
For HDF5/JLD2/Zarr backends, use [`write_operator`](@ref).

For [`W90Dat`](@ref):
- `write_w90_tb(file, tbdat, W90Dat())` writes only `prefix_tb.dat`.
- `write_w90_tb(file, tbdat, wsvec, W90Dat())` writes both `prefix_tb.dat` and `prefix_wsvec.dat`.
"""
function write_w90_tb end

function write_w90_tb(file::AbstractString, tbdat::TbDat)
    write_w90_tb_dat(file, tbdat)
    return nothing
end

function write_w90_tb(file::AbstractString, tbdat::TbDat, wsvec::WsvecDat)
    write_w90_tb_dat(file, tbdat)
    write_w90_wsvec_dat(_wsvec_path_from_tb(file), wsvec)
    return nothing
end

function _wsvec_path_from_tb(tbpath::AbstractString)
    if endswith(tbpath, "_tb.dat")
        return tbpath[1:(end - 7)] * "_wsvec.dat"
    else
        error("Not a valid tb.dat filename: $tbpath")
    end
end

function write_w90_tb(file::AbstractString, pack::OperatorPack)
    tbdat = TbDat(pack)
    wsvec = WsvecDat(pack)
    write_w90_tb(file, tbdat, wsvec)
    return nothing
end

function RvectorReducer(tbdat::TbDat)
    # Just the WsRvectorReducer
    return RvectorReducer(tbdat.Rvectors, tbdat.Rdegens)
end

function RvectorReducer(tbdat::TbDat, wsvec::WsvecDat)
    tbdat.Rvectors == wsvec.Rvectors ||
        error("Rvectors in tb.dat and wsvec.dat are inconsistent")

    if wsvec.mdrs
        return RvectorReducer(tbdat.Rvectors, tbdat.Rdegens, wsvec.Tvectors, wsvec.Tdegens)
    end
    return RvectorReducer(tbdat.Rvectors, tbdat.Rdegens)
end

function pack(tbdat::TbDat, wsvec::WsvecDat)
    reducer = RvectorReducer(tbdat, wsvec)
    return pack(tbdat, reducer)
end

function pack(tbdat::TbDat)
    reducer = RvectorReducer(tbdat)
    return pack(tbdat, reducer)
end

function pack(tbdat::TbDat, reducer::AbstractRvectorReducer)
    Rvectors = Vec3{Int}.(reducer.Rvectors)

    ops = OrderedDict(
        map([:H, :rx, :ry, :rz]) do name
            string(name) => reducer(getfield(tbdat, name))
        end,
    )

    Tr = eltype(tbdat.lattice)
    return OperatorPack(tbdat.header, Mat3{Tr}(tbdat.lattice), Rvectors, ops)
end

function TbDat(pack::OperatorPack)
    Tr = eltype(pack.lattice)
    Tc = Complex{Tr}
    Rvectors = Vec3{Int}.(pack.Rvectors)
    Rdegens = ones(Int, pack.n_Rvecs)

    H, rx, ry, rz = map(["H", "rx", "ry", "rz"]) do name
        haskey(pack.operators, name) || error("missing operator `$name`")
        op = pack.operators[name]
        if eltype(eltype(op)) <: Real
            [Tc.(O) for O in op]
        else
            op
        end
    end

    return TbDat(pack.header, Mat3{Tr}(pack.lattice), Rvectors, Rdegens, H, rx, ry, rz)
end

function WsvecDat(pack::OperatorPack)
    # Reduced OperatorPack does not preserve MDRS mapping, so export as non-MDRS.
    return WsvecDat(pack.header, Vec3{Int}.(pack.Rvectors), pack.n_wann)
end
