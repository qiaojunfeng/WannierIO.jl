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

    # length-`n_kpts` vector
    kpoints::Vector{Vec3{T}}

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
    # length-`n_kpts` vector, each element is a length-`n_bands` vector of bool
    dis_bands::Vector{BitVector}

    # 1D int array, size: n_kpts
    # n_dimfrozen:: Vector{Int}

    # u_matrix_opt, length-`n_kpts` vector, each elment has size: num_bands x num_wann
    Uᵈ::Vector{Matrix{Complex{T}}}

    # u_matrix, length-`n_kpts` vector, each element has size: num_wann x num_wann
    U::Vector{Matrix{Complex{T}}}

    # m_matrix, length-`n_kpts` vector of length-`n_bvecs` vector, each element is
    # matrix with size: num_wann x num_wann
    M::Vector{Vector{Matrix{Complex{T}}}}

    # wannier_centres, length-`num_wann` vector
    r::Vector{Vec3{T}}

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
    lattice::Mat3,
    recip_lattice::Mat3,
    kgrid::Vec3{Int},
    kpoints::Vector,
    checkpoint::String,
    have_disentangled::Bool,
    ΩI::Real,
    dis_bands::Vector{BitVector},
    Uᵈ::Vector,
    U::Vector,
    M::Vector,
    r::Vector,
    ω::Vector,
)
    if have_disentangled
        n_bands = size(Uᵈ[1], 1)
    else
        n_bands = size(U[1], 1)
    end

    n_exclude_bands = length(exclude_bands)
    n_kpts = length(M)
    n_bvecs = length(M[1])
    n_wann = size(U[1], 1)

    if have_disentangled
        n_dis = zeros(Int, n_kpts)
        for ik in 1:n_kpts
            n_dis[ik] = count(dis_bands[ik])
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
        collect.(Uᵈ),
        collect.(U),
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

    kpoints = zeros(Vec3{Float64}, n_kpts)
    for ik in 1:n_kpts
        kpoints[ik] = Vec3(parse.(Float64, split(srline()))...)
    end

    n_bvecs = parse(Int, srline())

    n_wann = parse(Int, srline())

    checkpoint = String(srline())

    # 1 -> True, 0 -> False
    have_disentangled = Bool(parse(Int, srline()))

    if have_disentangled
        # omega_invariant
        ΩI = parse(Float64, srline())

        dis_bands = [falses(n_bands) for _ in 1:n_kpts]
        for ik in 1:n_kpts
            for ib in 1:n_bands
                # 1 -> True, 0 -> False
                dis_bands[ik][ib] = Bool(parse(Int, srline()))
            end
        end

        n_dis = zeros(Int, n_kpts)
        for ik in 1:n_kpts
            n_dis[ik] = parse(Int, srline())
            @assert n_dis[ik] == count(dis_bands[ik])
        end

        # u_matrix_opt
        Uᵈ = [zeros(ComplexF64, n_bands, n_wann) for _ in 1:n_kpts]
        for ik in 1:n_kpts
            for iw in 1:n_wann
                for ib in 1:n_bands
                    vals = parse.(Float64, split(srline()))
                    Uᵈ[ik][ib, iw] = vals[1] + im * vals[2]
                end
            end
        end

    else
        ΩI = -1.0
        dis_bands = BitVector[]
        n_dis = Int[]
        Uᵈ = Matrix{ComplexF64}[]
    end

    # u_matrix
    U = [zeros(ComplexF64, n_wann, n_wann) for _ in 1:n_kpts]
    for ik in 1:n_kpts
        for iw in 1:n_wann
            for ib in 1:n_wann
                vals = parse.(Float64, split(srline()))
                U[ik][ib, iw] = vals[1] + im * vals[2]
            end
        end
    end

    #  m_matrix
    M = [[zeros(ComplexF64, n_wann, n_wann) for _ in 1:n_bvecs] for _ in 1:n_kpts]
    for ik in 1:n_kpts
        for inn in 1:n_bvecs
            for iw in 1:n_wann
                for ib in 1:n_wann
                    vals = parse.(Float64, split(srline()))
                    M[ik][inn][ib, iw] = vals[1] + im * vals[2]
                end
            end
        end
    end

    # wannier_centres
    r = zeros(Vec3{Float64}, n_wann)
    for iw in 1:n_wann
        r[iw] = Vec3(parse.(Float64, split(srline()))...)
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
        dis_bands = [tmp[:, i] for i in 1:n_kpts]
        n_dis = zeros(Int, n_kpts)
        n_dis .= read(io, (Tint, n_kpts))
        for ik in 1:n_kpts
            @assert n_dis[ik] == count(dis_bands[ik])
        end

        # u_matrix_opt
        U_tmp = zeros(ComplexF64, n_bands, n_wann, n_kpts)
        read(io, U_tmp)
        Uᵈ = [U_tmp[:, :, ik] for ik in 1:n_kpts]

    else
        ΩI = -1.0
        dis_bands = BitVector[]
        n_dis = Int[]
        Uᵈ = Matrix{ComplexF64}[]
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
        [U[:, :, ik] for ik in 1:n_kpts],
        [[M[:, :, ib, ik] for ib in 1:n_bvecs] for ik in 1:n_kpts],
        [Vec3(r[:, iw]) for iw in 1:n_wann],
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
                @printf(io, "%d\n", chk.dis_bands[ik][ib])
            end
        end

        for ik in 1:n_kpts
            @printf(io, "%d\n", chk.n_dis[ik])
        end

        # u_matrix_opt
        for ik in 1:n_kpts
            for iw in 1:n_wann
                for ib in 1:n_bands
                    v = chk.Uᵈ[ik][ib, iw]
                    @printf(io, "%25.17f %25.17f\n", real(v), imag(v))
                end
            end
        end
    end

    # u_matrix
    for ik in 1:n_kpts
        for iw in 1:n_wann
            for ib in 1:n_wann
                v = chk.U[ik][ib, iw]
                @printf(io, "%25.17f %25.17f\n", real(v), imag(v))
            end
        end
    end

    #  m_matrix
    for ik in 1:n_kpts
        for inn in 1:n_bvecs
            for iw in 1:n_wann
                for ib in 1:n_wann
                    v = chk.M[ik][inn][ib, iw]
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
    # size = 3 * n_kpts
    kpoints = hcat(chk.kpoints...)
    write(io, Float64.(kpoints))

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
        write(io, Tint.(hcat(chk.dis_bands)))

        write(io, Tint.(chk.n_dis))

        # u_matrix_opt
        write(io, ComplexF64.(cat(chk.Uᵈ...; dims=3)))
    end

    # u_matrix
    write(io, ComplexF64.(cat(chk.U...; dims=3)))

    #  m_matrix
    M = zeros(ComplexF64, n_wann, n_wann, n_bvecs, n_kpts)
    for ik in 1:n_kpts
        for inn in 1:n_bvecs
            M[:, :, inn, ik] = chk.M[ik][inn]
        end
    end
    write(io, M)

    # wannier_centres
    r = hcat(chk.r...)
    write(io, Float64.(r))

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
    get_U(chk)

Extract `U` matrices from `Chk`.
"""
function get_U(chk::Chk)
    if !chk.have_disentangled
        # Return deepcopy for safety, so that chk.U is not modified
        return deepcopy(chk.U)
    end

    Uᵈ = get_Udis(chk)

    return map(zip(Uᵈ, chk.U)) do (U1, U2)
        U1 * U2
    end
end

"""
    get_Udis(chk)

Extract `U` matrices for disentanglement from `Chk`.
"""
function get_Udis(chk::Chk)
    n_kpts = chk.n_kpts
    n_bands = chk.n_bands
    n_wann = chk.n_wann

    T = eltype(chk.U[1])
    if !chk.have_disentangled
        return [diagm(0 => ones(T, n_wann)) for _ in 1:n_kpts]
    end

    # Need to permute wavefunctions since Uᵈ is stored in a way that
    # the bands taking part in disentanglement are in the first few rows.
    # Construct identity matrix
    Iᵏ = Matrix{T}(I, n_bands, n_bands)

    return map(1:n_kpts) do ik
        # sortperm is stable, and
        # need descending order (dis bands at the front)
        p = sortperm(chk.dis_bands[ik]; order=Base.Order.Reverse)
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
        Iᵏ[:, p] * chk.Uᵈ[ik]
    end
end

"""
Compare two `Chk` objects.

Used in tests.
"""
function Base.isapprox(a::Chk, b::Chk)
    for f in propertynames(a)
        va = getfield(a, f)
        vb = getfield(b, f)

        if va isa String
            va == vb || return false
        elseif va isa Vector
            if eltype(va) isa BitVector
                all(va .== vb) || return false
            else
                all(va .≈ vb) || return false
            end
        else
            va ≈ vb || return false
        end
    end
    return true
end
