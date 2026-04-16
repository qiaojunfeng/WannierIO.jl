export read_chk, write_chk, gauge_matrices, gauge_matrices_dis

"""
Struct for storing matrices in `prefix.chk` file.

$(TYPEDEF)

One-to-one mapping to the wannier90 `chk` file, but renaming the variable
names so that they are consistent with the rest of the code.

# Fields

$(FIELDS)
"""
struct Chk{T <: Real}
    "The header line, usually contains date and time"
    header::String

    """number of bands, can be auto set in constructor according to dimensions
    of other variables"""
    n_bands::Int

    "number of excluded bands, can be auto set in constructor"
    n_exclude_bands::Int

    """Indices of excluded bands, starts from 1.
    Vector of integers, size: `n_exclude_bands`"""
    exclude_bands::Vector{Int}

    "Matrix of size 3 x 3, each column is a lattice vector in Å unit"
    lattice::Mat3{T}

    "Matrix of size 3 x 3, each column is a reciprocal lattice vector in Å⁻¹ unit"
    recip_lattice::Mat3{T}

    "number of kpoints, can be auto set in constructor"
    n_kpts::Int

    "dimensions of kpoint grid, 3 integers"
    kgrid::Vec3{Int}

    "kpoint coordinates, fractional, length-`n_kpts` vector"
    kpoints::Vector{Vec3{T}}

    "number of b-vectors, can be auto set in constructor"
    n_bvecs::Int

    "number of Wannier functions, can be auto set in constructor"
    n_wann::Int

    """a string to indicate the current step (after disentanglement, after
    maximal localization, ...) in wannier90"""
    checkpoint::String

    """Have finished disentanglement or not"""
    have_disentangled::Bool

    "Omega invariant part of MV spreads, in Å² unit"
    ΩI::T

    """Indices of bands taking part in disentanglement, not frozen bands!
    stored as a dense `n_bands × n_kpts` boolean matrix.

    This is needed since W90 puts all the disentanglement bands
    in the first several rows of `Udis`,
    (and the first few columns of `Udis` are the frozen bands)
    so directly multiplying eigenvalues e.g.
    `(Udis * U)' * diag(eigenvalues) * (Udis * U)` is wrong!
    """
    dis_bands::BitMatrix

    """number of bands taking part in disentanglement at each kpoint.
    can be auto set in constructor from `dis_bands`"""
    n_dis::Vector{Int}

    """Semi-unitary matrix for disentanglement with size `n_bands × n_wann × n_kpts`,
    i.e., the `u_matrix_opt` in wannier90."""
    Udis::Array{Complex{T}, 3}

    """Unitary matrix for maximal localization with size `n_wann × n_wann × n_kpts`,
    i.e., the `u_matrix` in wannier90.
    The abbreviation `ml` stands for maximal localization, so as to
    differentiate from the (combined) unitary matrix `U = Udis * Uml`."""
    Uml::Array{Complex{T}, 3}

    """Wannier-gauge overlap matrix with size `n_wann × n_wann × n_bvecs × n_kpts`,
    i.e., the `m_matrix` in wannier90"""
    M::Array{Complex{T}, 4}

    """Wannier function centers, length-`n_wann` vector, Cartesian coordinates
    in Å unit, i.e., the `wannier_centres` variable in wannier90"""
    r::Vector{Vec3{T}}

    """Wannier function spreads, length-`n_wann` vector, Å² unit,
    i.e., the `wannier_spreads` variable in wannier90"""
    ω::Vector{T}
end

# I still save these size in the Chk struct, as these number are indeed written
# in the chk file, the Chk struct tries to minic that layout.
# These functions are provided to have a unified API.
n_bands(chk::Chk) = chk.n_bands
n_wannier(chk::Chk) = chk.n_wann
n_bvectors(chk::Chk) = chk.n_bvecs
n_kpoints(chk::Chk) = chk.n_kpts

function Base.show(io::IO, chk::Chk)
    return print(io, "Chk(n_kpts=$(chk.n_kpts), n_bands=$(chk.n_bands), n_wann=$(chk.n_wann))")
end

function Base.show(io::IO, ::MIME"text/plain", chk::Chk)
    return print(
        io,
        """Chk(
          header: $(chk.header)
          checkpoint: $(chk.checkpoint)
          n_kpts: $(chk.n_kpts)
          n_bands: $(chk.n_bands)
          n_wann: $(chk.n_wann)
          n_exclude_bands: $(chk.n_exclude_bands)
          have_disentangled: $(chk.have_disentangled)
          ΩI: $(chk.ΩI)
        )""",
    )
end

"""
    $(SIGNATURES)

Convenience constructor of [`Chk`](@ref) struct that auto set some fields.
"""
function Chk(
        header::AbstractString,
        exclude_bands::AbstractVector{Int},
        lattice::AbstractMatrix,
        recip_lattice::AbstractMatrix,
        kgrid::AbstractVector{<:Integer},
        kpoints::AbstractVector,
        checkpoint::AbstractString,
        have_disentangled::Bool,
        ΩI::Real,
        dis_bands::AbstractMatrix{Bool},
        Udis::AbstractArray,
        Uml::AbstractArray,
        M::AbstractArray,
        r::AbstractVector,
        ω::AbstractVector,
    )
    dis_bands = BitMatrix(dis_bands)
    Udis = Array(Udis)
    Uml = Array(Uml)
    M = Array(M)

    n_bands = have_disentangled ? size(Udis, 1) : size(Uml, 1)
    n_exclude_bands = length(exclude_bands)
    n_wann, n_wann2, n_bvecs, n_kpts = size(M)
    n_wann == n_wann2 || error("M is not square in Wannier indices")
    size(Uml, 1) == size(Uml, 2) || error("Uml is not square")
    n_wann == size(Uml, 1) || error("inconsistent n_wann between M and Uml")

    if have_disentangled
        n_dis = zeros(Int, n_kpts)
        for ik in 1:n_kpts
            n_dis[ik] = count(@view dis_bands[:, ik])
        end
    else
        n_dis = zeros(Int, 0)
    end

    return Chk(
        header,
        n_bands,
        n_exclude_bands,
        collect(exclude_bands),
        mat3(lattice),
        mat3(recip_lattice),
        n_kpts,
        kgrid,
        kpoints,
        n_bvecs,
        n_wann,
        checkpoint,
        have_disentangled,
        ΩI,
        dis_bands,
        n_dis,
        Udis,
        Uml,
        M,
        r,
        ω,
    )
end

"""
    read_chk(filename)
    read_chk(file, ::FortranText)
    read_chk(file, ::FortranBinary)

Read wannier90 `chk` checkpoint file.

# Arguments
- `file`: The name of the input file, or an `IO`.

Similar to [`read_amn`](@ref), the 1st version auto detect `chk` file format
(binary or text) and read it.
"""
function read_chk(io::IO, ::FortranText)
    header = String(readstrip(io))

    n_bands = parse(Int, readstrip(io))

    n_exclude_bands = parse(Int, readstrip(io))

    exclude_bands = zeros(Int, n_exclude_bands)
    for i in 1:n_exclude_bands
        exclude_bands[i] = parse(Int, readstrip(io))
    end

    # Each column is a lattice vector
    # but W90 writes x components first, then y, z. NOT a1 first, then a2, a3.
    line = parse_vector(readstrip(io))
    lattice = Mat3{Float64}(reshape(line, (3, 3))')

    # Each column is a lattice vector
    line = parse_vector(readstrip(io))
    recip_lattice = Mat3{Float64}(reshape(line, (3, 3))')

    n_kpts = parse(Int, readstrip(io))

    kgrid = Vec3{Int}(parse.(Int, split(readstrip(io))))

    kpoints = zeros(Vec3{Float64}, n_kpts)
    for ik in 1:n_kpts
        kpoints[ik] = parse_vector(readstrip(io))
    end

    n_bvecs = parse(Int, readstrip(io))
    n_wann = parse(Int, readstrip(io))
    checkpoint = String(readstrip(io))

    # 1 -> True, 0 -> False
    have_disentangled = Bool(parse(Int, readstrip(io)))

    if have_disentangled
        # omega_invariant
        ΩI = parse(Float64, readstrip(io))

        dis_bands = falses(n_bands, n_kpts)
        for ik in 1:n_kpts
            for ib in 1:n_bands
                # 1 -> True, 0 -> False
                dis_bands[ib, ik] = Bool(parse(Int, readstrip(io)))
            end
        end

        n_dis = zeros(Int, n_kpts)
        for ik in 1:n_kpts
            n_dis[ik] = parse(Int, readstrip(io))
            n_dis[ik] == count(@view dis_bands[:, ik]) ||
                error("Inconsistent number of disentangled bands")
        end

        # u_matrix_opt
        Udis = zeros(ComplexF64, n_bands, n_wann, n_kpts)
        for ik in 1:n_kpts
            for iw in 1:n_wann
                for ib in 1:n_bands
                    vals = parse_vector(readstrip(io))
                    Udis[ib, iw, ik] = vals[1] + im * vals[2]
                end
            end
        end
    else
        ΩI = -1.0
        dis_bands = falses(0, 0)
        n_dis = Int[]
        Udis = zeros(ComplexF64, 0, 0, 0)
    end

    # u_matrix
    Uml = zeros(ComplexF64, n_wann, n_wann, n_kpts)
    for ik in 1:n_kpts
        for iw in 1:n_wann
            for ib in 1:n_wann
                vals = parse_vector(readstrip(io))
                Uml[ib, iw, ik] = vals[1] + im * vals[2]
            end
        end
    end

    #  m_matrix
    M = zeros(ComplexF64, n_wann, n_wann, n_bvecs, n_kpts)
    for ik in 1:n_kpts
        for inn in 1:n_bvecs
            for iw in 1:n_wann
                for ib in 1:n_wann
                    vals = parse_vector(readstrip(io))
                    M[ib, iw, inn, ik] = vals[1] + im * vals[2]
                end
            end
        end
    end

    # wannier_centres
    r = zeros(Vec3{Float64}, n_wann)
    for iw in 1:n_wann
        r[iw] = parse_vector(readstrip(io))
    end

    # wannier_spreads
    ω = zeros(Float64, n_wann)
    for iw in 1:n_wann
        ω[iw] = parse(Float64, readstrip(io))
    end

    return Chk(
        header,
        exclude_bands,
        lattice,
        recip_lattice,
        kgrid,
        kpoints,
        checkpoint,
        have_disentangled,
        ΩI,
        dis_bands,
        Udis,
        Uml,
        M,
        r,
        ω,
    )
end

function read_chk(filename::AbstractString, ::FortranText)
    return open(filename) do io
        read_chk(io, FortranText())
    end
end

function read_chk(io::FortranFile, ::FortranBinary)

    # strip and read line
    header_len = 33
    header = trimstring(read(io, FString{header_len}))

    # gfortran default integer is 4 bytes
    Tint = Int32
    n_bands = read(io, Tint)

    n_exclude_bands = read(io, Tint)

    exclude_bands = zeros(Int, n_exclude_bands)
    exclude_bands .= read(io, (Tint, n_exclude_bands))
    # I don't know why but I need to skip a record here
    # probably because the FortranFiles.read does not handle 0-length arrays correctly
    n_exclude_bands == 0 && read(io)

    # Each column is a lattice vector
    # but W90 writes x components first, then y, z. Not a1 first, then a2, a3.
    lattice = read(io, (Float64, 3, 3))
    lattice = Mat3{Float64}(lattice')

    # Each column is a lattice vector
    recip_lattice = read(io, (Float64, 3, 3))
    recip_lattice = Mat3{Float64}(recip_lattice')

    n_kpts = read(io, Tint)

    kgrid = Vec3{Int}(read(io, (Tint, 3)))

    kpoints = zeros(Vec3{Float64}, n_kpts)
    read(io, kpoints)

    n_bvecs = read(io, Tint)

    n_wann = read(io, Tint)

    checkpoint = trimstring(read(io, FString{20}))

    # treat Bool as Int32
    # 1 -> true, 0 -> false
    have_disentangled = parse_bool(read(io, Tint))

    if have_disentangled
        # omega_invariant
        ΩI = read(io, Float64)

        tmp = parse_bool.(read(io, (Tint, n_bands, n_kpts)))
        dis_bands = BitMatrix(tmp)
        n_dis = zeros(Int, n_kpts)
        n_dis .= read(io, (Tint, n_kpts))
        for ik in 1:n_kpts
            n_dis[ik] == count(@view dis_bands[:, ik]) ||
                error("Inconsistent number of disentangled bands")
        end

        # u_matrix_opt
        Udis = zeros(ComplexF64, n_bands, n_wann, n_kpts)
        read(io, Udis)

    else
        ΩI = -1.0
        dis_bands = falses(0, 0)
        n_dis = Int[]
        Udis = zeros(ComplexF64, 0, 0, 0)
    end

    # u_matrix
    Uml = zeros(ComplexF64, n_wann, n_wann, n_kpts)
    read(io, Uml)

    #  m_matrix
    M = zeros(ComplexF64, n_wann, n_wann, n_bvecs, n_kpts)
    read(io, M)

    # wannier_centres
    r = zeros(Float64, 3, n_wann)
    read(io, r)

    # wannier_spreads
    ω = zeros(Float64, n_wann)
    read(io, ω)

    close(io)

    return Chk(
        header,
        exclude_bands,
        lattice,
        recip_lattice,
        kgrid,
        kpoints,
        checkpoint,
        have_disentangled,
        ΩI,
        dis_bands,
        Udis,
        Uml,
        M,
        [vec3(r[:, iw]) for iw in 1:n_wann],
        ω,
    )
end

function read_chk(filename::AbstractString, ::FortranBinary)
    io = FortranFile(filename)
    return read_chk(io, FortranBinary())
end

function read_chk(filename::AbstractString)
    format = detect_fortran_format(filename)
    return read_chk(filename, format)
end

"""
    write_chk(filename, chk::Chk; binary=false)
    write_chk(file, chk::Chk, ::FortranText)
    write_chk(file, chk::Chk, ::FortranBinary)

Write wannier90 `chk` file.

Similar to [`write_amn`](@ref), the 1st version is a convenience wrapper.

# Arguments
- `file`: The name of the output file, or an `IO`.
- `chk`: the [`Chk`](@ref) struct

# Keyword arguments
- `binary`: write as Fortran binary file or not. Although wannier90 default
    is Fortran binary format, here the default is `false` since Fortran binary
    depends on compiler and platform, so it is not guaranteed to always work.
"""
function write_chk end

function write_chk(io::IO, chk::Chk, ::FortranText)
    n_bands = chk.n_bands
    n_wann = chk.n_wann
    n_kpts = chk.n_kpts
    n_bvecs = chk.n_bvecs

    # Write formatted chk file
    @printf(io, "%33s\n", chk.header)

    @printf(io, "%d\n", n_bands)

    @printf(io, "%d\n", chk.n_exclude_bands)

    if chk.n_exclude_bands > 0
        for i in 1:(chk.n_exclude_bands)
            @printf(io, "%d\n", chk.exclude_bands[i])
        end
    end

    # Each column is a lattice vector
    # but W90 writes x components first, then y, z. Not a1 first, then a2, a3.
    for v in reshape(chk.lattice', 9)
        @printf(io, "%25.17f", v)
    end
    @printf(io, "\n")

    # Each column is a lattice vector
    for v in reshape(chk.recip_lattice', 9)
        @printf(io, "%25.17f", v)
    end
    @printf(io, "\n")

    @printf(io, "%d\n", n_kpts)

    @printf(io, "%d %d %d\n", chk.kgrid...)

    for kpt in chk.kpoints
        @printf(io, "%25.17f %25.17f %25.17f\n", kpt...)
    end

    @printf(io, "%d\n", n_bvecs)

    @printf(io, "%d\n", n_wann)

    # left-justified
    @printf(io, "%-20s\n", chk.checkpoint)

    # 1 -> True, 0 -> False
    # v = chk.have_disentangled ? 1 : 0
    @printf(io, "%d\n", chk.have_disentangled)

    if chk.have_disentangled
        # omega_invariant
        @printf(io, "%25.17f\n", chk.ΩI)

        for ik in 1:n_kpts
            for ib in 1:n_bands
                # 1 -> True, 0 -> False
                @printf(io, "%d\n", chk.dis_bands[ib, ik])
            end
        end

        for ik in 1:n_kpts
            @printf(io, "%d\n", chk.n_dis[ik])
        end

        # u_matrix_opt
        for ik in 1:n_kpts
            for iw in 1:n_wann
                for ib in 1:n_bands
                    v = chk.Udis[ib, iw, ik]
                    @printf(io, "%25.17f %25.17f\n", real(v), imag(v))
                end
            end
        end
    end

    # u_matrix
    for ik in 1:n_kpts
        for iw in 1:n_wann
            for ib in 1:n_wann
                v = chk.Uml[ib, iw, ik]
                @printf(io, "%25.17f %25.17f\n", real(v), imag(v))
            end
        end
    end

    #  m_matrix
    for ik in 1:n_kpts
        for inn in 1:n_bvecs
            for iw in 1:n_wann
                for ib in 1:n_wann
                    v = chk.M[ib, iw, inn, ik]
                    @printf(io, "%25.17f %25.17f\n", real(v), imag(v))
                end
            end
        end
    end

    # wannier_centres
    for iw in 1:n_wann
        @printf(io, "%25.17f %25.17f %25.17f\n", chk.r[iw]...)
    end

    # wannier_spreads
    for iw in 1:n_wann
        @printf(io, "%25.17f\n", chk.ω[iw])
    end
    return nothing
end

function write_chk(filename::AbstractString, chk::Chk, ::FortranText)
    open(filename, "w") do io
        write_chk(io, chk, FortranText())
    end
    return nothing
end

function write_chk(io::FortranFile, chk::Chk, ::FortranBinary)
    n_bands = chk.n_bands
    n_wann = chk.n_wann
    n_kpts = chk.n_kpts
    n_bvecs = chk.n_bvecs

    write(io, FString(33, chk.header))

    # gfortran default integer is 4 bytes
    Tint = Int32
    write(io, Tint(n_bands))

    write(io, Tint(chk.n_exclude_bands))

    write(io, Vector{Tint}(chk.exclude_bands))

    # Each column is a lattice vector
    # but W90 writes x components first, then y, z. Not a1 first, then a2, a3.
    write(io, Vector{Float64}(reshape(chk.lattice', 9)))

    # Each column is a lattice vector
    write(io, Vector{Float64}(reshape(chk.recip_lattice', 9)))

    write(io, Tint(n_kpts))

    write(io, Vector{Tint}(chk.kgrid))

    # size = 3 * n_kpts
    # Use `reduce(hcat, kpoints)`, which is much faster than `hcat(kpoints...)`
    kpoints = reduce(hcat, chk.kpoints)
    write(io, Matrix{Float64}(kpoints))

    write(io, Tint(n_bvecs))

    write(io, Tint(n_wann))

    # left-justified
    write(io, FString(20, chk.checkpoint))

    # true -> 1, false -> 0
    write(io, Tint(chk.have_disentangled))

    if chk.have_disentangled
        # omega_invariant
        write(io, Float64(chk.ΩI))

        # true -> 1, false -> 0
        write(io, Matrix{Tint}(chk.dis_bands))

        write(io, Vector{Tint}(chk.n_dis))

        # u_matrix_opt
        write(io, chk.Udis)
    end

    # u_matrix
    write(io, chk.Uml)

    #  m_matrix
    write(io, chk.M)

    # wannier_centres
    r = reduce(hcat, chk.r)
    write(io, Matrix{Float64}(r))

    # wannier_spreads
    write(io, Vector{Float64}(chk.ω))

    close(io)
    return nothing
end

function write_chk(filename::AbstractString, chk::Chk, ::FortranBinary)
    io = FortranFile(filename, "w")
    return write_chk(io, chk, FortranBinary())
end

function write_chk(filename::AbstractString, chk::Chk; binary = false)
    format = fortran_format(; binary)
    return write_chk(filename, chk, format)
end

"""
    $(SIGNATURES)

Extract the combined `U = Udis * Uml` matrices from `Chk`.
"""
function gauge_matrices(chk::Chk)
    Uml = gauge_matrices_ml(chk)
    if !chk.have_disentangled
        return Uml
    end

    Udis = gauge_matrices_dis(chk)
    U = zeros(eltype(Uml), chk.n_bands, chk.n_wann, chk.n_kpts)
    for ik in 1:chk.n_kpts
        U[:, :, ik] = Udis[:, :, ik] * Uml[:, :, ik]
    end
    return U
end

function gauge_matrices_ml(chk::Chk)
    # Return deepcopy for safety, so that chk.Uml is not modified
    return deepcopy(chk.Uml)
end

"""
    $(SIGNATURES)

Extract disentanglement `Udis` matrices from `Chk`.
"""
function gauge_matrices_dis(chk::Chk)
    n_kpts = chk.n_kpts
    n_bands = chk.n_bands
    n_wann = chk.n_wann

    T = eltype(chk.Uml)
    if !chk.have_disentangled
        U = zeros(T, n_bands, n_wann, n_kpts)
        for ik in 1:n_kpts
            for iw in 1:n_wann
                U[iw, iw, ik] = 1
            end
        end
        return U
    end

    # Need to permute wavefunctions since Udis is stored in a way that
    # the bands taking part in disentanglement are in the first few rows.
    # Construct identity matrix
    Iᵏ = Matrix{T}(I, n_bands, n_bands)
    U = zeros(T, n_bands, n_wann, n_kpts)
    for ik in 1:n_kpts
        # sortperm is stable, and
        # need descending order (dis bands at the front)
        p = sortperm(view(chk.dis_bands, :, ik); order = Base.Order.Reverse)
        # usually we don't need this permutation, but if
        # 1. the dis_win_min > minimum(E), then these below
        #    dis_win_min bands are shifted to the last rows of Udis
        # 2. use projectability disentanglement, then
        #    there might be cases that the lower (bonding) and
        #    higher (anti-bonding) bands participate in disentanglement,
        #    but some low-projectability bands are excluded from
        #    disentanglement, then these low-proj bands are shifted to
        #    the last rows of Udis
        # so we need to permute the Bloch states before multiplying Udis
        # chk.Udis: semi-unitary matrices from disentanglement
        # chk.Uml: unitary matrices from maximal localization
        U[:, :, ik] = Iᵏ[:, p] * chk.Udis[:, :, ik]
    end
    return U
end

"""
    $(SIGNATURES)

Compare two `Chk` objects.
"""
function Base.isapprox(a::Chk, b::Chk)
    return _isapprox(a, b)
end
