@testitem "WsRvectorReducer" begin
    using LazyArtifacts

    tbdat = WannierIO.read_w90_tb_dat(artifact"Si2_valence/outputs/WS/Si2_valence_tb.dat")

    reducer = WannierIO.WsRvectorReducer(tbdat.Rvectors, tbdat.Rdegens)
    # Test the exported function `RvectorReducer` works
    reducer2 = RvectorReducer(tbdat.Rvectors, tbdat.Rdegens)
    @test reducer == reducer2

    H1 = reducer(tbdat.H)
    rx1 = reducer(tbdat.rx)
    ry1 = reducer(tbdat.ry)
    rz1 = reducer(tbdat.rz)

    ref = WannierIO.read_w90_tb_dat(
        artifact"Si2_valence/outputs/WS/reduced_Rvectors/Si2_valence_tb.dat"
    )

    atol = 1.0e-10
    @test reducer.Rvectors == ref.Rvectors
    @test reducer.degens == tbdat.Rdegens
    @test isapprox(H1, ref.H; atol)
    @test isapprox(rx1, ref.rx; atol)
    @test isapprox(ry1, ref.ry; atol)
    @test isapprox(rz1, ref.rz; atol)
end

@testitem "MdrsRvectorReducer" begin
    using LazyArtifacts

    tbdat = WannierIO.read_w90_tb_dat(artifact"Si2_valence/outputs/MDRS/Si2_valence_tb.dat")
    wsvec = WannierIO.read_w90_wsvec_dat(
        artifact"Si2_valence/outputs/MDRS/Si2_valence_wsvec.dat"
    )

    reducer = WannierIO.MdrsRvectorReducer(
        tbdat.Rvectors, tbdat.Rdegens, wsvec.Tvectors, wsvec.Tdegens
    )
    # Test the exported function `RvectorReducer` works
    reducer2 = RvectorReducer(tbdat.Rvectors, tbdat.Rdegens, wsvec.Tvectors, wsvec.Tdegens)
    @test reducer == reducer2

    H1 = reducer(tbdat.H)
    rx1 = reducer(tbdat.rx)
    ry1 = reducer(tbdat.ry)
    rz1 = reducer(tbdat.rz)

    ref_tbdat = WannierIO.read_w90_tb_dat(
        artifact"Si2_valence/outputs/MDRS/reduced_Rvectors/Si2_valence_tb.dat"
    )
    ref_wsvec = WannierIO.read_w90_wsvec_dat(
        artifact"Si2_valence/outputs/MDRS/Si2_valence_wsvec.dat"
    )

    # The ref tb.dat was generated with a different ordering of the new R vectors
    idx = map(reducer.Rvectors) do R
        findfirst(==(R), ref_tbdat.Rvectors)
    end
    ref_Rvectors = ref_tbdat.Rvectors[idx]
    ref_H = ref_tbdat.H[:, :, idx]
    ref_rx = ref_tbdat.rx[:, :, idx]
    ref_ry = ref_tbdat.ry[:, :, idx]
    ref_rz = ref_tbdat.rz[:, :, idx]

    atol = 1.0e-10
    @test reducer.Rvectors == ref_Rvectors
    @test isapprox(H1, ref_H; atol)
    @test isapprox(rx1, ref_rx; atol)
    @test isapprox(ry1, ref_ry; atol)
    @test isapprox(rz1, ref_rz; atol)
end
