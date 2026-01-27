@testitem "read isym" begin
    using LazyArtifacts
    sym = read_isym(artifact"Si2_hse/Si2.isym")

    @test sym.n_symops == 96
    @test sym.spinors == false

    @test sym.symops[end].R == [0 -1 0; 1 1 1; 0 0 -1]
    @test sym.symops[end].t ≈ [-1 / 4, -1 / 4, -1 / 4]
    @test sym.symops[end].time_reversal == true
    @test sym.symops[end].u == zeros(ComplexF64, 2, 2)
    @test sym.symops[end].isym == 96
    @test sym.symops[end].isym_inv == 91

    @test sym.nkpts_ibz == 29
    @test sym.kpoints_ibz[end] ≈ [1 / 4, -1 / 2, -1 / 4]

    @test sym.n_bands == 16
    @test sym.n_repmat_band == 340

    @test sym.repmat_band[end].ik_ibz == 29
    @test sym.repmat_band[end].isym == 82
    @test sym.repmat_band[end].d[1, 2] ≈ -0.024533833915732 - 0.093336077514977im

    @test sym.n_wann == 8

    @test sym.repmat_wann[end].isym == 96
    @test sym.repmat_wann[end].D[1, 5] ≈ 0.999999999999999
end

@testitem "build_mapping_ik_isym" begin
    using LazyArtifacts
    sym = read_isym(artifact"Si2_hse/Si2.isym")
    mapping = WannierIO.build_mapping_ik_isym(sym.repmat_band; sym.nkpts_ibz, sym.n_symops)

    @test mapping[1][1] == 1
    @test mapping[29][82] == length(sym.repmat_band)
end
