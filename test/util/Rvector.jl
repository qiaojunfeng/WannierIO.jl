@testitem "WsRvectorReducer" begin
    using LazyArtifacts

    tbdat = WannierIO.read_w90_tb_dat(artifact"Si2_valence/outputs/WS/Si2_valence_tb.dat")

    reducer = WannierIO.WsRvectorReducer(tbdat.Rvectors, tbdat.Rdegens)
    H1 = reducer(tbdat.H)
    r_x1 = reducer(tbdat.r_x)
    r_y1 = reducer(tbdat.r_y)
    r_z1 = reducer(tbdat.r_z)

    ref = WannierIO.read_w90_tb_dat(
        artifact"Si2_valence/outputs/WS/reduced_Rvectors/Si2_valence_tb.dat"
    )

    atol = 1e-10
    @test reducer.Rvectors == ref.Rvectors
    @test reducer.degens == tbdat.Rdegens
    @test isapprox(H1, ref.H; atol)
    @test isapprox(r_x1, ref.r_x; atol)
    @test isapprox(r_y1, ref.r_y; atol)
    @test isapprox(r_z1, ref.r_z; atol)
end

@testitem "MdrsRvectorReducer" begin
    using LazyArtifacts

    tbdat = WannierIO.read_w90_tb_dat(artifact"Si2_valence/outputs/MDRS/Si2_valence_tb.dat")
    wsvec = WannierIO.read_w90_wsvec(
        artifact"Si2_valence/outputs/MDRS/Si2_valence_wsvec.dat"
    )

    reducer = WannierIO.MdrsRvectorReducer(
        tbdat.Rvectors, tbdat.Rdegens, wsvec.Tvectors, wsvec.Tdegens
    )
    H1 = reducer(tbdat.H)
    r_x1 = reducer(tbdat.r_x)
    r_y1 = reducer(tbdat.r_y)
    r_z1 = reducer(tbdat.r_z)

    ref_tbdat = WannierIO.read_w90_tb_dat(
        artifact"Si2_valence/outputs/MDRS/reduced_Rvectors/Si2_valence_tb.dat"
    )
    ref_wsvec = WannierIO.read_w90_wsvec(
        artifact"Si2_valence/outputs/MDRS/Si2_valence_wsvec.dat"
    )

    # The ref tb.dat was generated with a different ordering of the new R vectors
    idx = map(reducer.Rvectors) do R
        findfirst(==(R), ref_tbdat.Rvectors)
    end
    ref_Rvectors = ref_tbdat.Rvectors[idx]
    ref_H = ref_tbdat.H[idx]
    ref_r_x = ref_tbdat.r_x[idx]
    ref_r_y = ref_tbdat.r_y[idx]
    ref_r_z = ref_tbdat.r_z[idx]

    atol = 1e-10
    @test reducer.Rvectors == ref_Rvectors
    @test isapprox(H1, ref_H; atol)
    @test isapprox(r_x1, ref_r_x; atol)
    @test isapprox(r_y1, ref_r_y; atol)
    @test isapprox(r_z1, ref_r_z; atol)
end
