export read_w90_tb, write_w90_tb

"""
    read_w90_tb(file, format)
    read_w90_tb(file)

Read Wannier90 tight-binding data.

The Rvectors and degeneracies will be reduced by [`RvectorReducer`](@ref)
automatically, resulting in a [`OperatorPack`](@ref) with just the Rvectors
and the operator matrices defined on them.

If a corresponding `wsvec.dat` file is found, the MDRS Tvectors and degeneracies
will also be reduced automatically.

As a consequence, the Wannier interpolation will be just a simple inverse
Fourier transform without any complications from R-degeneracies or MDRS
T vectors and degeneracies.

See also [`read_w90_tb_dat`](@ref) and [`read_w90_wsvec`](@ref) for reading just
the `prefix_tb.dat` and `prefix_wsvec.dat` files without any reduction.

The format is one of:

- [`W90Dat`](@ref): native Wannier90 `.dat` text file
- [`HDF5Format`](@ref): HDF5 (requires loading `HDF5`)
- [`JLD2Format`](@ref): JLD2 (requires loading `JLD2`)
- [`ZarrFormat`](@ref): Zarr (requires loading `Zarr`)

When called without a format argument the format is inferred from the file
extension via [`detect_w90dat_format`](@ref).
"""
function read_w90_tb end

function read_w90_tb(file, ::W90Dat)
    tbdat = read_w90_tb_dat(file)
    wsvec_path = _wsvec_path_from_tb(file)
    if isfile(wsvec_path)
        wsvec = read_w90_wsvec(wsvec_path)
        return pack(tbdat, wsvec)
    else
        @warn "No corresponding wsvec file found for tb file at $file. \
               Reading tb.dat without wsvec information."
    end
    return pack(tbdat)
end

"""
    read_w90_tb backend fallback — override by loading HDF5/JLD2/Zarr.
"""
function read_w90_tb(file, fmt::AbstractFileFormat)
    error("Format `$(format_name(fmt))` requires loading the corresponding package. \
           See the WannierIO documentation for supported formats.")
end

function read_w90_tb(file)
    return read_w90_tb(file, detect_w90dat_format(file))
end

"""
    write_w90_tb(file, tbdat, format; kwargs...)
    write_w90_tb(file, tbdat)

    write_w90_tb(file, tbdat, wsvec, format; kwargs...)
    write_w90_tb(file, tbdat, wsvec)

    write_w90_tb(file, pack, format; kwargs...)
    write_w90_tb(file, pack)

Write Wannier90 tight-binding data. The format is one of:

- [`W90Dat`](@ref): native Wannier90 `.dat` text file
- [`HDF5Format`](@ref): HDF5 — keywords: `atol`, `index_type`, `value_type`, `deflate`
- [`JLD2Format`](@ref): JLD2 — keywords: `atol`, `index_type`, `value_type`, `compress`
- [`ZarrFormat`](@ref): Zarr — keywords: `atol`, `index_type`, `value_type`, `clevel`

When called without a format argument, the format is inferred from the file
extension via [`detect_w90dat_format`](@ref).

The `kwargs` are passed to the corresponding backend and control the sparsification
and compression of the TB data. See the documentation of each format for details.

For [`W90Dat`](@ref):
- `write_w90_tb(file, tbdat, W90Dat())` writes only `prefix_tb.dat`.
- `write_w90_tb(file, tbdat, wsvec, W90Dat())` writes both `prefix_tb.dat` and `prefix_wsvec.dat`.
- `write_w90_tb(file, pack::OperatorPack, W90Dat())` writes one `prefix_tb.dat`
    with unit R-degeneracies.
"""
function write_w90_tb end

function write_w90_tb(file::AbstractString, tbdat::TbDat, ::W90Dat)
    write_w90_tb_dat(file, tbdat)
    return nothing
end

function write_w90_tb(file::AbstractString, tbdat::TbDat, wsvec::WsvecDat, ::W90Dat)
    write_w90_tb_dat(file, tbdat)
    write_w90_wsvec(_wsvec_path_from_tb(file), wsvec)
    return nothing
end

function write_w90_tb(file::AbstractString, pack::OperatorPack, ::W90Dat; kwargs...)
    write_w90_tb_dat(file, TbDat(pack))
    return nothing
end

"""
    write_w90_tb backend fallback — override by loading HDF5/JLD2/Zarr.
"""
function write_w90_tb(
    file::AbstractString, pack::OperatorPack, fmt::AbstractFileFormat; kwargs...
)
    error("Format `$(format_name(fmt))` requires loading the corresponding package. \
           See the `WannierIO` documentation for supported formats.")
end

# Convenience: TbDat input is converted to OperatorPack first.
function write_w90_tb(
    file::AbstractString, tbdat::TbDat, fmt::AbstractFileFormat; kwargs...
)
    return write_w90_tb(file, pack(tbdat), fmt; kwargs...)
end

# Convenience: WsvecDat + TbDat input applies WS/MDRS reduction first.
function write_w90_tb(
    file::AbstractString, tbdat::TbDat, wsvec::WsvecDat, fmt::AbstractFileFormat; kwargs...
)
    return write_w90_tb(file, pack(tbdat, wsvec), fmt; kwargs...)
end

function write_w90_tb(file::AbstractString, tbdat::TbDat; kwargs...)
    return write_w90_tb(file, tbdat, detect_w90dat_format(file); kwargs...)
end

function write_w90_tb(file::AbstractString, tbdat::TbDat, wsvec::WsvecDat; kwargs...)
    return write_w90_tb(file, tbdat, wsvec, detect_w90dat_format(file); kwargs...)
end

function write_w90_tb(file::AbstractString, pack::OperatorPack; kwargs...)
    return write_w90_tb(file, pack, detect_w90dat_format(file); kwargs...)
end

function _wsvec_path_from_tb(tbpath::AbstractString)
    if endswith(tbpath, "_tb.dat")
        return tbpath[1:(end - 7)] * "_wsvec.dat"
    else
        error("Not a valid tb.dat filename: $tbpath")
    end
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
        map([:H, :r_x, :r_y, :r_z]) do name
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

    H, r_x, r_y, r_z = map(["H", "r_x", "r_y", "r_z"]) do name
        haskey(pack.operators, name) || error("missing operator `$name`")
        op = pack.operators[name]
        if eltype(eltype(op)) <: Real
            [Tc.(O) for O in op]
        else
            op
        end
    end

    return TbDat(pack.header, Mat3{Tr}(pack.lattice), Rvectors, Rdegens, H, r_x, r_y, r_z)
end

"""
    $(SIGNATURES)

Convert dense `TbDat` operators into a [`SparseOperatorPack`](@ref).
"""
function sparsify(tbdat::TbDat, opt::SparseOption=SparseOption())
    sparsify(pack(tbdat), opt)
end

function sparsify(tbdat::TbDat, wsvec::WsvecDat, opt::SparseOption=SparseOption())
    sparsify(pack(tbdat, wsvec), opt)
end
