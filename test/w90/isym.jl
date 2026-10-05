@testitem "read isym raw" begin
    using LazyArtifacts, LinearAlgebra
    sym = read_isym_raw(artifact"Si2_hse/Si2.isym")

    @test sym.n_symops == 96
    @test sym.spinors == false

    @test sym.symops[end].s == [0 1 0; -1 1 0; 0 1 -1]
    @test sym.symops[end].ft ≈ [-1 / 4, -1 / 4, -1 / 4]
    @test sym.symops[end].t_rev == true
    @test sym.symops[end].u == Matrix{ComplexF64}(I, 2, 2)
    @test sym.symops[end].isym == 96
    @test sym.symops[end].invs == 91

    @test sym.nkpts_ibz == 29
    @test sym.kpoints_ibz[end] ≈ [1 / 4, -1 / 2, -1 / 4]

    @test sym.n_bands == 16
    @test length(sym.repmat_band) == 340

    @test sym.repmat_band[end].ik_ibz == 29
    @test sym.repmat_band[end].ig == 82
    @test sym.repmat_band[end].d[1, 2] ≈ -0.024533833915732 - 0.093336077514977im

    @test sym.n_wann == 8

    @test length(sym.repmat_wann) == 96
    @test sym.repmat_wann[end].D[1, 5] ≈ 0.999999999999999
end

@testitem "read spinor isym matrix" begin
    io = IOBuffer(
        """
        spinor parser regression
        1 1
        identity
        1 0 0
        0 1 0
        0 0 1
        0.0 0.0 0.0
        0
        1.0 2.0
        3.0 4.0
        5.0 6.0
        7.0 8.0
        1

        K points
        1
        0.0 0.0 0.0

        Representation matrix of G_k
        1 1
        1 1 1
        1 1 1.0 0.0

        Rotation matrix of Wannier functions
        1
        1 1
        1 1 1.0 0.0
        """
    )

    sym = read_isym_raw(io)
    @test sym.spinors
    # pw2wannier90 writes the four SU(2) entries in row-major order.
    @test sym.symops[1].u == ComplexF64[
        1 + 2im 3 + 4im
        5 + 6im 7 + 8im
    ]
    @test sym.symops[1].invs == 1
end

@testitem "standardize isym" begin
    using LazyArtifacts, LinearAlgebra
    using WannierIO: Mat3
    raw = read_isym_raw(artifact"Si2_hse/Si2.isym")
    sym = standardize(raw)

    @test sym.n_symops == raw.n_symops
    @test sym.spinors == raw.spinors
    @test sym.nkpts_ibz == raw.nkpts_ibz
    @test sym.kpoints_ibz == raw.kpoints_ibz
    @test sym.n_bands == raw.n_bands
    @test sym.n_wann == raw.n_wann

    # the little-group representations are unchanged
    @test sym.littlegroup_reps == raw.repmat_band

    for (op, rawop) in zip(sym.symops, raw.symops)
        # k-space rotation is the file's s matrix
        @test op.Wk == rawop.s
        # W = transpose(inv(s)), integer
        @test op.W * transpose(rawop.s) == Mat3{Int}(I)
        # v = W * ft
        @test op.v ≈ op.W * rawop.ft
        @test op.time_reversal == rawop.t_rev
        @test op.ig_inv == rawop.invs
    end

    # orbital_reps[ig] stores D(g_ig): the raw file stores D(g_ig^{-1})
    # at index ig, so the standardized entry is its inverse (the adjoint
    # for a unitary operation, the transpose for an antiunitary one)
    @test any(op -> op.time_reversal, sym.symops)
    @test all(op.ig == ig for (ig, op) in enumerate(sym.symops))
    for ig in 1:sym.n_symops
        Draw = raw.repmat_wann[ig].D
        @test sym.orbital_reps[ig].D == (raw.symops[ig].t_rev ? transpose(Draw) : Draw')
    end

    # read_isym is the standardized read
    @test read_isym(artifact"Si2_hse/Si2.isym").symops[end].W == sym.symops[end].W
end

@testitem "SymOp one-line display" begin
    using LazyArtifacts
    symops = read_isym(artifact"Si2_hse/Si2.isym").symops

    @test repr(symops[1]) == "SymOp(1: identity)"
    @test repr(symops[5]) == "SymOp(5: 180 deg rotation - cart. axis [1,1,0], v = [0.25, 0.25, -0.75])"
    # the -0.0 components of a pure rotation do not print as a translation
    @test repr(symops[2]) == "SymOp(2: 180 deg rotation - cart. axis [0,0,1])"
    @test repr(symops[49]) == "SymOp(49: identity, +T)"
    # one line per operation inside a container
    @test countlines(IOBuffer(sprint(show, MIME"text/plain"(), symops))) == 1 + length(symops)

    bare = WannierIO.SymOp("", symops[1].W, symops[1].v, symops[1].Wk, false, symops[1].u, 1, 1)
    @test repr(bare) == "SymOp(1)"
end

@testitem "orbital representations are stored at the index the file gives" begin
    # two operations whose orbital entries the file lists in reverse order
    isym_text(entries) = """
    order regression
    2 0
    identity
    1 0 0
    0 1 0
    0 0 1
    0.0 0.0 0.0
    0
    1
    inversion
    -1 0 0
    0 -1 0
    0 0 -1
    0.0 0.0 0.0
    0
    2

    K points
    1
    0.0 0.0 0.0

    Representation matrix of G_k
    1 2
    1 1 1
    1 1 1.0 0.0
    1 2 1
    1 1 -1.0 0.0

    Rotation matrix of Wannier functions
    1
    $entries
    """
    raw = read_isym_raw(IOBuffer(isym_text("2 1\n1 1 -1.0 0.0\n1 1\n1 1 1.0 0.0")))
    @test raw.repmat_wann[1].D == fill(1.0 + 0im, 1, 1)
    @test raw.repmat_wann[2].D == fill(-1.0 + 0im, 1, 1)
    @test raw.repmat_band[2].ig == 2

    @test_throws ErrorException read_isym_raw(IOBuffer(isym_text("1 1\n1 1 1.0 0.0\n1 1\n1 1 1.0 0.0")))
    @test_throws ErrorException read_isym_raw(IOBuffer(isym_text("3 1\n1 1 1.0 0.0\n1 1\n1 1 1.0 0.0")))
end
