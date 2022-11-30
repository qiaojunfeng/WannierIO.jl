using Printf: @printf

export read_chk, write_chk, get_U, get_Udis

"""
Struct for storing matrices in `seedname.chk` file.
"""
struct Chk{T<:Real}
    header::String

    # These commented variables are written inside chk file,
    # but I move them to the end of the struct and use a constructor
    # to set them according to the shape of the corresponding arrays.

    # n_bands:: Int
    # n_exclude_bands:: Int

    # 1D array, size: num_exclude_bands
    exclude_bands::Vector{Int}

    # 2D array, size: 3 x 3, each column is a lattice vector
    lattice::Mat3{T}

    # 2D array, size: 3 x 3, each column is a lattice vector
    recip_lattice::Mat3{T}

    # n_kpts: Int

    # List of 3 int
    kgrid::Vec3{Int}

    # 2D array, size: 3 x n_kpts, each column is a vector
    kpoints::Matrix{T}

    # n_bvecs:: Int

    # n_wann: int

    checkpoint::String

    have_disentangled::Bool

    # omega invariant
    ΩI::T

    # Bands taking part in disentanglement, not frozen bands!
    # This is needed since W90 puts all the disentanglement bands
    # in the first several rows of Uᵈ,
    # (and the first few columns of Uᵈ are the frozen bands)
    # so directly multiplying eigenvalues e.g.
    # (Uᵈ * U)' * diag(eigenvalues) * (Uᵈ * U) is wrong!
    # 2D bool array, size: n_bands x n_kpts
    dis_bands::BitMatrix

    # 1D int array, size: n_kpts
    # n_dimfrozen:: Vector{Int}

    # u_matrix_opt, 3D array, size: num_bands x num_wann x num_kpts
    Uᵈ::Array{Complex{T},3}

    # u_matrix, 3D array, size: num_wann x num_wann x num_kpts
    U::Array{Complex{T},3}

    # m_matrix, 4D array, size:num_wann x num_wann x n_bvecs x num_kpts
    M::Array{Complex{T},4}

    # wannier_centres, 2D array, size: 3 x num_wann
    r::Matrix{T}

    # wannier_spreads, 1D array, size: num_wann
    ω::Vector{T}

    # these variables are auto-set in constructor
    n_bands::Int
    n_exclude_bands::Int
    n_kpts::Int
    n_bvecs::Int
    n_wann::Int
    n_dis::Vector{Int}
end

function Chk(
    header::String,
    exclude_bands::Vector{Int},
    lattice::Mat3{T},
    recip_lattice::Mat3{T},
    kgrid::Vec3{Int},
    kpoints::Matrix{T},
    checkpoint::String,
    have_disentangled::Bool,
    ΩI::T,
    dis_bands::BitMatrix,
    Uᵈ::Array{Complex{T},3},
    U::Array{Complex{T},3},
    M::Array{Complex{T},4},
    r::Matrix{T},
    ω::Vector{T},
) where {T<:Real}
    if have_disentangled
        n_bands = size(Uᵈ, 1)
    else
        n_bands = size(U, 1)
    end

    n_exclude_bands = length(exclude_bands)

    n_kpts = size(M, 4)

    n_bvecs = size(M, 3)

    n_wann = size(U, 1)

    if have_disentangled
        n_dis = zeros(Int, n_kpts)
        for ik in 1:n_kpts
            n_dis[ik] = count(dis_bands[:, ik])
        end
    else
        n_dis = zeros(Int, 0)
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
        Uᵈ,
        U,
        M,
        r,
        ω,
        n_bands,
        n_exclude_bands,
        n_kpts,
        n_bvecs,
        n_wann,
        n_dis,
    )
end

"""
Read formatted `seedname.chk` file.
"""
function _read_chk_fmt(filename::AbstractString)
    io = open(filename)

    # strip and read line
    srline() = strip(readline(io))

    # Read formatted chk file
    header = String(srline())

    n_bands = parse(Int, srline())

    n_exclude_bands = parse(Int, srline())

    exclude_bands = zeros(Int, n_exclude_bands)

    if n_exclude_bands > 0
        for i in 1:n_exclude_bands
            exclude_bands[i] = parse(Int, srline())
        end
    end

    # Each column is a lattice vector
    # but W90 writes x components first, then y, z. Not a1 first, then a2, a3.
    line = parse.(Float64, split(srline()))
    lattice = Mat3{Float64}(reshape(line, (3, 3))')

    # Each column is a lattice vector
    line = parse.(Float64, split(srline()))
    recip_lattice = Mat3{Float64}(reshape(line, (3, 3))')

    n_kpts = parse(Int, srline())

    kgrid = Vec3{Int}(parse.(Int, split(srline())))

    kpoints = zeros(Float64, 3, n_kpts)
    for ik in 1:n_kpts
        kpoints[:, ik] = parse.(Float64, split(srline()))
    end

    n_bvecs = parse(Int, srline())

    n_wann = parse(Int, srline())

    checkpoint = String(srline())

    # 1 -> True, 0 -> False
    have_disentangled = Bool(parse(Int, srline()))

    if have_disentangled
        # omega_invariant
        ΩI = parse(Float64, srline())

        dis_bands = falses(n_bands, n_kpts)
        for ik in 1:n_kpts
            for ib in 1:n_bands
                # 1 -> True, 0 -> False
                dis_bands[ib, ik] = Bool(parse(Int, srline()))
            end
        end

        n_dis = zeros(Int, n_kpts)
        for ik in 1:n_kpts
            n_dis[ik] = parse(Int, srline())
            @assert n_dis[ik] == count(dis_bands[:, ik])
        end

        # u_matrix_opt
        Uᵈ = zeros(ComplexF64, n_bands, n_wann, n_kpts)
        for ik in 1:n_kpts
            for iw in 1:n_wann
                for ib in 1:n_bands
                    vals = parse.(Float64, split(srline()))
                    Uᵈ[ib, iw, ik] = vals[1] + im * vals[2]
                end
            end
        end

    else
        ΩI = -1.0
        dis_bands = falses(0, 0)
        n_dis = zeros(Int, 0)
        Uᵈ = zeros(ComplexF64, 0, 0, 0)
    end

    # u_matrix
    U = zeros(ComplexF64, n_wann, n_wann, n_kpts)
    for ik in 1:n_kpts
        for iw in 1:n_wann
            for ib in 1:n_wann
                vals = parse.(Float64, split(srline()))
                U[ib, iw, ik] = vals[1] + im * vals[2]
            end
        end
    end

    #  m_matrix
    M = zeros(ComplexF64, n_wann, n_wann, n_bvecs, n_kpts)
    for ik in 1:n_kpts
        for inn in 1:n_bvecs
            for iw in 1:n_wann
                for ib in 1:n_wann
                    vals = parse.(Float64, split(srline()))
                    M[ib, iw, inn, ik] = vals[1] + im * vals[2]
                end
            end
        end
    end

    # wannier_centres
    r = zeros(Float64, 3, n_wann)
    for iw in 1:n_wann
        r[:, iw] = parse.(Float64, split(srline()))
    end

    # wannier_spreads
    ω = zeros(Float64, n_wann)
    for iw in 1:n_wann
        ω[iw] = parse(Float64, srline())
    end

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
        Uᵈ,
        U,
        M,
        r,
        ω,
    )
end

"""
Read unformatted `seedname.chk` file.
"""
function _read_chk_bin(filename::AbstractString)
    io = FortranFile(filename)

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

    kpoints = zeros(Float64, 3, n_kpts)
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

        dis_bands = falses(n_bands, n_kpts)
        dis_bands .= parse_bool.(read(io, (Tint, n_bands, n_kpts)))

        n_dis = zeros(Int, n_kpts)
        n_dis .= read(io, (Tint, n_kpts))
        for ik in 1:n_kpts
            @assert n_dis[ik] == count(dis_bands[:, ik])
        end

        # u_matrix_opt
        Uᵈ = zeros(ComplexF64, n_bands, n_wann, n_kpts)
        read(io, Uᵈ)
    else
        ΩI = -1.0
        dis_bands = falses(0, 0)
        n_dis = zeros(Int, 0)
        Uᵈ = zeros(ComplexF64, 0, 0, 0)
    end

    # u_matrix
    U = zeros(ComplexF64, n_wann, n_wann, n_kpts)
    read(io, U)

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
        Uᵈ,
        U,
        M,
        r,
        ω,
    )
end

"""
    read_chk(filename::AbstractString)

Read formatted or binary `chk` file.
"""
function read_chk(filename::AbstractString)
    @info "Reading chk file:" filename
    println()

    if isbinary(filename)
        return _read_chk_bin(filename)
    else
        return _read_chk_fmt(filename)
    end
end

"""
Write formatted `chk` file.
"""
function _write_chk_fmt(filename::AbstractString, chk::Chk)
    io = open(filename, "w")

    n_bands = chk.n_bands
    n_wann = chk.n_wann
    n_kpts = chk.n_kpts
    n_bvecs = chk.n_bvecs

    # Read formatted chk file
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

    for ik in 1:n_kpts
        @printf(io, "%25.17f %25.17f %25.17f\n", chk.kpoints[:, ik]...)
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
                    v = chk.Uᵈ[ib, iw, ik]
                    @printf(io, "%25.17f %25.17f\n", real(v), imag(v))
                end
            end
        end
    end

    # u_matrix
    for ik in 1:n_kpts
        for iw in 1:n_wann
            for ib in 1:n_wann
                v = chk.U[ib, iw, ik]
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
        @printf(io, "%25.17f %25.17f %25.17f\n", chk.r[:, iw]...)
    end

    # wannier_spreads
    for iw in 1:n_wann
        @printf(io, "%25.17f\n", chk.ω[iw])
    end

    return close(io)
end

"""
Write unformatted `chk` file.
"""
function _write_chk_bin(filename::AbstractString, chk::Chk)
    io = FortranFile(filename, "w")

    n_bands = chk.n_bands
    n_wann = chk.n_wann
    n_kpts = chk.n_kpts
    n_bvecs = chk.n_bvecs

    write(io, FString(33, chk.header))

    # gfortran default integer is 4 bytes
    Tint = Int32
    write(io, Tint(n_bands))

    write(io, Tint(chk.n_exclude_bands))

    write(io, Tint.(chk.exclude_bands))

    # Each column is a lattice vector
    # but W90 writes x components first, then y, z. Not a1 first, then a2, a3.
    write(io, Float64.(reshape(chk.lattice', 9)))

    # Each column is a lattice vector
    write(io, Float64.(reshape(chk.recip_lattice', 9)))

    write(io, Tint(n_kpts))

    write(io, Tint.(chk.kgrid))

    write(io, Float64.(chk.kpoints))

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
        write(io, Tint.(chk.dis_bands))

        write(io, Tint.(chk.n_dis))

        # u_matrix_opt
        write(io, ComplexF64.(chk.Uᵈ))
    end

    # u_matrix
    write(io, ComplexF64.(chk.U))

    #  m_matrix
    write(io, ComplexF64.(chk.M))

    # wannier_centres
    write(io, Float64.(chk.r))

    # wannier_spreads
    write(io, Float64.(chk.ω))

    return close(io)
end

"""
    write_chk(filename, chk::Chk)

Write formatted or binary `chk` file.

# Keyword arguments
- `binary`: write as Fortran unformatted file
"""
function write_chk(filename::AbstractString, chk::Chk; binary::Bool=false)
    if binary
        _write_chk_bin(filename, chk)
    else
        _write_chk_fmt(filename, chk)
    end

    @info "Written to file: $(filename)"
    println()
    return nothing
end

"""
    get_A(chk::Chk)

Extract `U` matrices from `Chk`.
"""
function get_U(chk::Chk)
    n_kpts = chk.n_kpts
    n_bands = chk.n_bands
    n_wann = chk.n_wann

    U = similar(chk.U, n_bands, n_wann, n_kpts)

    if !chk.have_disentangled
        U .= chk.U
        return U
    end

    Uᵈ = get_Udis(chk)

    for ik in 1:n_kpts
        # Uᵈ: semi-unitary matrices from disentanglement
        # chk.U: unitary matrices from maximal localization
        U[:, :, ik] = Uᵈ[:, :, ik] * chk.U[:, :, ik]
    end

    return U
end

"""
    get_Udis(chk::Chk)

Extract `A` matrices for disentanglement from `Chk`.
"""
function get_Udis(chk::Chk)
    n_kpts = chk.n_kpts
    n_bands = chk.n_bands
    n_wann = chk.n_wann

    if !chk.have_disentangled
        return eyes_A(eltype(U), n_bands, n_kpts)
    end

    U = similar(chk.U, n_bands, n_wann, n_kpts)

    # need to permute wavefunctions since Uᵈ is stored in a way that
    # the bands taking part in disentanglement are in the first few rows
    Iᵏ = Matrix{eltype(U)}(I, n_bands, n_bands)

    for ik in 1:n_kpts
        # sortperm is stable, and
        # need descending order (dis bands at the front)
        p = sortperm(chk.dis_bands[:, ik]; order=Base.Order.Reverse)
        # usually we don't need this permutation, but if
        # 1. the dis_win_min > minimum(E), then these below
        #    dis_win_min bands are shifted to the last rows of Uᵈ
        # 2. use projectability disentanglement, then
        #    there might be cases that the lower (bonding) and
        #    higher (anti-bonding) bands participate in disentanglement,
        #    but some low-projectability bands are excluded from
        #    disentanglement, then these low-proj bands are shifted to
        #    the last rows of Uᵈ
        # so we need to permute the Bloch states before multiplying Uᵈ
        # chk.Uᵈ: semi-unitary matrices from disentanglement
        # chk.U: unitary matrices from maximal localization
        U[:, :, ik] = Iᵏ[:, p] * chk.Uᵈ[:, :, ik]
    end

    return U
end
