using Printf: @printf

export read_chk, write_chk, get_U, get_Udis

"""
Struct for storing matrices in `prefix.chk` file.

$(TYPEDEF)

One-to-one mapping to the wannier90 `chk` file, but renaming the variable
names so that they are consistent with the rest of the code.

# Fields

$(FIELDS)
"""
struct Chk{T<:Real}
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

    """
    Indices of bands taking part in disentanglement, not frozen bands!
    length-`n_kpts` vector, each element is a length-`n_bands` vector of bool.

    This is needed since W90 puts all the disentanglement bands
    in the first several rows of `Udis`,
    (and the first few columns of `Udis` are the frozen bands)
    so directly multiplying eigenvalues e.g.
    `(Udis * U)' * diag(eigenvalues) * (Udis * U)` is wrong!
    """
    dis_bands::Vector{BitVector}

    """number of bands taking part in disentanglement at each kpoint.
    can be auto set in constructor from `dis_bands`"""
    n_dis::Vector{Int}

    """Semi-unitary matrix for disentanglement,
    length-`n_kpts` vector, each elment has size: `n_bands` x `n_wann`,
    i.e., the `u_matrix_opt` in wannier90"""
    Udis::Vector{Matrix{Complex{T}}}

    """Unitary matrix for maximal localization,
    length-`n_kpts` vector, each element has size: `n_wann` x `n_wann`,
    i.e., the `u_matrix` in wannier90.
    The abbreviation `ml` stands for maximal localization, so as to
    differentiate from the (combined) unitary matrix `U = Udis * Uml`."""
    Uml::Vector{Matrix{Complex{T}}}

    """Wannier-gauge overlap matrix,
    length-`n_kpts` vector of length-`n_bvecs` vector, each element is
    a matrix of size `n_wann` x `n_wann`,
    i.e., the `m_matrix` in wannier90"""
    M::Vector{Vector{Matrix{Complex{T}}}}

    """Wannier function centers, length-`n_wann` vector, Cartesian coordinates
    in Å unit, i.e., the `wannier_centres` variable in wannier90"""
    r::Vector{Vec3{T}}

    """Wannier function spreads, length-`n_wann` vector, Å² unit,
    i.e., the `wannier_spreads` variable in wannier90"""
    ω::Vector{T}
end

"""
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
    dis_bands::AbstractVector{BitVector},
    Udis::AbstractVector,
    Uml::AbstractVector,
    M::AbstractVector,
    r::AbstractVector,
    ω::AbstractVector,
)
    if have_disentangled
        @assert length(Udis) > 0 "empty Udis"
        n_bands = size(Udis[1], 1)
    else
        @assert length(Uml) > 0 "empty Uml"
        n_bands = size(Uml[1], 1)
    end

    n_exclude_bands = length(exclude_bands)
    n_kpts = length(M)
    @assert n_kpts > 0 "empty M"
    n_bvecs = length(M[1])
    @assert length(Uml) > 0 "empty Uml"
    n_wann = size(Uml[1], 1)

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
        n_bands,
        n_exclude_bands,
        exclude_bands,
        Mat3(lattice),
        Mat3(recip_lattice),
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
        collect.(Udis),
        collect.(Uml),
        M,
        r,
        ω,
    )
end

"""
    read_chk(filename)
    read_chk(filename, ::FortranText)
    read_chk(filename, ::FortranBinary)

Read wannier90 `chk` checkpoint file.

Similar to [`read_amn`](@ref), the 1st version auto detect `chk` file format
(binary or text) and read it.
"""
function read_chk(filename::AbstractString, ::FortranText)
    chk = open(filename) do io
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
        # but W90 writes x components first, then y, z. NOT a1 first, then a2, a3.
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
            Udis = [zeros(ComplexF64, n_bands, n_wann) for _ in 1:n_kpts]
            for ik in 1:n_kpts
                for iw in 1:n_wann
                    for ib in 1:n_bands
                        vals = parse.(Float64, split(srline()))
                        Udis[ik][ib, iw] = vals[1] + im * vals[2]
                    end
                end
            end

        else
            ΩI = -1.0
            dis_bands = BitVector[]
            n_dis = Int[]
            Udis = Matrix{ComplexF64}[]
        end

        # u_matrix
        Uml = [zeros(ComplexF64, n_wann, n_wann) for _ in 1:n_kpts]
        for ik in 1:n_kpts
            for iw in 1:n_wann
                for ib in 1:n_wann
                    vals = parse.(Float64, split(srline()))
                    Uml[ik][ib, iw] = vals[1] + im * vals[2]
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
    return chk
end

function read_chk(filename::AbstractString, ::FortranBinary)
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
        Udis = [U_tmp[:, :, ik] for ik in 1:n_kpts]

    else
        ΩI = -1.0
        dis_bands = BitVector[]
        n_dis = Int[]
        Udis = Matrix{ComplexF64}[]
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
        [Uml[:, :, ik] for ik in 1:n_kpts],
        [[M[:, :, ib, ik] for ib in 1:n_bvecs] for ik in 1:n_kpts],
        [Vec3(r[:, iw]) for iw in 1:n_wann],
        ω,
    )
end

function read_chk(filename::AbstractString)
    if isbinary(filename)
        format = FortranBinary()
    else
        format = FortranText()
    end
    chk = read_chk(filename, format)

    n_kpts = chk.n_kpts
    n_bands = chk.n_bands
    n_wann = chk.n_wann
    @info "Reading chk file" filename n_kpts n_bands n_wann
    return chk
end

"""
    write_chk(filename, chk::Chk; binary=false)
    write_chk(filename, chk::Chk, ::FortranText)
    write_chk(filename, chk::Chk, ::FortranBinary)

Write wannier90 `chk` file.

Similar to [`write_amn`](@ref), the 1st version is a convenience wrapper.

# Keyword arguments
- `binary`: write as Fortran binary file or not. Although wannier90 default
    is Fortran binary format, here the default is `false` since Fortran binary
    depends on compiler and platform, so it is not guaranteed to always work.
"""
function write_chk end

function write_chk(filename::AbstractString, chk::Chk, ::FortranText)
    open(filename, "w") do io
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
                        v = chk.Udis[ik][ib, iw]
                        @printf(io, "%25.17f %25.17f\n", real(v), imag(v))
                    end
                end
            end
        end

        # u_matrix
        for ik in 1:n_kpts
            for iw in 1:n_wann
                for ib in 1:n_wann
                    v = chk.Uml[ik][ib, iw]
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
    end
end

function write_chk(filename::AbstractString, chk::Chk, ::FortranBinary)
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

    # concatenate along dims=3
    cat3(A...) = cat(A...; dims=3)

    if chk.have_disentangled
        # omega_invariant
        write(io, Float64(chk.ΩI))

        # true -> 1, false -> 0
        write(io, Matrix{Tint}(reduce(hcat, chk.dis_bands)))

        write(io, Vector{Tint}(chk.n_dis))

        # u_matrix_opt
        write(io, Array{ComplexF64}(reduce(cat3, chk.Udis)))
    end

    # u_matrix
    write(io, Array{ComplexF64}(reduce(cat3, chk.Uml)))

    #  m_matrix
    M = zeros(ComplexF64, n_wann, n_wann, n_bvecs, n_kpts)
    for ik in 1:n_kpts
        for ibvec in 1:n_bvecs
            M[:, :, ibvec, ik] = chk.M[ik][ibvec]
        end
    end
    write(io, M)

    # wannier_centres
    r = reduce(hcat, chk.r)
    write(io, Matrix{Float64}(r))

    # wannier_spreads
    write(io, Vector{Float64}(chk.ω))

    close(io)
    return nothing
end

function write_chk(filename::AbstractString, chk::Chk; binary=false)
    n_kpts = chk.n_kpts
    n_bands = chk.n_bands
    n_wann = chk.n_wann
    @info "Writing chk file" filename n_kpts n_bands n_wann

    if binary
        format = FortranBinary()
    else
        format = FortranText()
    end
    return write_chk(filename, chk, format)
end

"""
    get_U(chk)

Extract `U` matrices from `Chk`.
"""
function get_U(chk::Chk)
    if !chk.have_disentangled
        # Return deepcopy for safety, so that chk.Uml is not modified
        return deepcopy(chk.Uml)
    end

    Udis = get_Udis(chk)
    return map(zip(Udis, chk.Uml)) do (d, m)
        d * m
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

    T = eltype(chk.Uml[1])
    if !chk.have_disentangled
        return [diagm(0 => ones(T, n_wann)) for _ in 1:n_kpts]
    end

    # Need to permute wavefunctions since Udis is stored in a way that
    # the bands taking part in disentanglement are in the first few rows.
    # Construct identity matrix
    Iᵏ = Matrix{T}(I, n_bands, n_bands)

    return map(1:n_kpts) do ik
        # sortperm is stable, and
        # need descending order (dis bands at the front)
        p = sortperm(chk.dis_bands[ik]; order=Base.Order.Reverse)
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
        Iᵏ[:, p] * chk.Udis[ik]
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
