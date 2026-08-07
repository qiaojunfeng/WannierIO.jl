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
    @test sym.repmat_band[end].isym == 82
    @test sym.repmat_band[end].d[1, 2] ≈ -0.024533833915732 - 0.093336077514977im

    @test sym.n_wann == 8

    @test sym.repmat_wann[end].isym == 96
    @test sym.repmat_wann[end].D[1, 5] ≈ 0.999999999999999
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
        @test op.isym_inv == rawop.invs
    end

    # orbital_reps[isym] stores D(g_isym): the raw file stores D(g_isym^{-1})
    # at index isym, so the standardized entry is the raw entry at invs(isym)
    for isym in 1:sym.n_symops
        @test sym.orbital_reps[isym].isym == isym
        @test sym.orbital_reps[isym].D == raw.repmat_wann[raw.symops[isym].invs].D
    end

    # read_isym is the standardized read
    @test read_isym(artifact"Si2_hse/Si2.isym").symops[end].W == sym.symops[end].W
end

@testitem "build_mapping_ik_isym" begin
    using LazyArtifacts
    sym = read_isym(artifact"Si2_hse/Si2.isym")
    mapping = WannierIO.build_mapping_ik_isym(
        sym.littlegroup_reps; sym.nkpts_ibz, sym.n_symops
    )

    @test mapping[1][1] == 1
    @test mapping[29][82] == length(sym.littlegroup_reps)
end
