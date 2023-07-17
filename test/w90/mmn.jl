
@testset "read/write mmn" begin
    M, kpb_k, kpb_G = read_mmn(artifact"Si2_valence/Si2_valence.mmn")

    @test length(M) == 216
    @test length(M[1]) == 8
    @test size(M[1][1]) == (4, 4)

    ref_M215_2 = ComplexF64[
        0.622355772788+0.779341618658im -0.015531265051+0.027176244708im -0.003508661654-0.015077859568im 0.021680417129+0.039444115733im
        -0.00654020394-0.043590860733im -0.716198599581+0.416420660955im -0.022335323703+0.073160414474im -0.011094081817+0.495937237686im
        0.000395190251-0.001890191419im -0.473340195142+0.082043373917im 0.346341288525-0.454328488051im 0.219116863121-0.552770061648im
        0.007094594145-0.002164553263im 0.15507683774+0.19419614694im 0.69754000164+0.104779280281im -0.363933244214+0.045670933434im
    ]

    @test M[215][2] ≈ ref_M215_2

    ref_kpb_k_3 = [2, 4, 9, 39, 46, 33, 183, 212]
    @test kpb_k[3] == ref_kpb_k_3

    ref_kpb_G_3 = WannierIO.Vec3[
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        [0, -1, 0],
        [-1, 0, 0],
        [-1, -1, 0],
    ]
    @test kpb_G[3] == ref_kpb_G_3

    tmpfile = tempname(; cleanup=true)
    write_mmn(tmpfile, M, kpb_k, kpb_G)
    M2, kpb_k2, kpb_G2 = read_mmn(tmpfile)

    @test M ≈ M2
    @test kpb_k ≈ kpb_k2
    @test kpb_G ≈ kpb_G2
end

@testset "read/write mmn binary" begin
    M, kpb_k, kpb_G = read_mmn(artifact"Si2_valence/Si2_valence.mmn")
    M1, kpb_k1, kpb_G1 = read_mmn(artifact"Si2_valence/reference/binary/Si2_valence.mmn")
    @test M ≈ M1
    @test kpb_k ≈ kpb_k1
    @test kpb_G ≈ kpb_G1

    tmpfile = tempname(; cleanup=true)
    write_mmn(tmpfile, M, kpb_k, kpb_G; binary=true)
    M2, kpb_k2, kpb_G2 = read_mmn(tmpfile)

    @test M ≈ M2
    @test kpb_k ≈ kpb_k2
    @test kpb_G ≈ kpb_G2
end
