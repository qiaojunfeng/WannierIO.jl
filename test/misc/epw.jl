@testitem "read EPW mmn" begin
    using LazyArtifacts
    M, kpb_k, kpb_G = read_mmn(artifact"BN/BN.mmn")
    n_kpts = length(M)
    n_bvecs = length(M[1])
    n_bands = size(M[1][1], 1)
    M1 = WannierIO.read_epw_mmn(artifact"BN/outputs/bn.mmn"; n_kpts, n_bvecs, n_bands)
    @test M ≈ M1
end

@testitem "read/write EPW ukk" begin
    using LazyArtifacts
    chk = read_chk(artifact"BN/outputs/bn.chk")
    alat = WannierIO.read_qe_xml(artifact"BN/outputs/bn.xml").alat
    ukk_ref = WannierIO.Ukk(chk, alat)

    ukk = WannierIO.read_epw_ukk(artifact"BN/outputs/bn.ukk")
    @test ukk ≈ ukk_ref

    tmpfile = tempname(; cleanup=true)
    WannierIO.write_epw_ukk(tmpfile, ukk)

    ukk2 = WannierIO.read_epw_ukk(tmpfile)
    @test ukk2 ≈ ukk
end
