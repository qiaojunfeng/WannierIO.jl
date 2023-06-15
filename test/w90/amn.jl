@testset "read/write amn" begin
    A = read_amn(artifact"Si2_valence/si2.amn")

    ref_A215 = ComplexF64[
        0.602605803842+0.466458612517im 0.709819979644+0.096505743176im 0.602605806825+0.466458612235im 0.271333571639+0.662974997524im
        -2.7652e-8+2.436e-9im -0.401551121567-0.311690904391im -1.7245e-8+1.8091e-8im -0.069156613881+0.503598920796im
        -0.009606395172-0.299196334599im 0.150309036867+0.2420558235im -0.009606388763-0.299196326647im -0.134471911066+0.251199304047im
        -0.392902559295+0.232501544615im 5.7803e-8-7.7606e-8im 0.39290241635-0.232501390525im 1.50756e-7-8.2322e-8im
    ]

    @test length(A) == 216
    @test size(A[1]) == (4, 4)
    @test A[215] ≈ ref_A215

    tmpfile = tempname(; cleanup=true)
    write_amn(tmpfile, A)

    A1 = read_amn(tmpfile)
    @test A ≈ A1
end

@testset "read/write amn binary" begin
    A = read_amn(artifact"Si2_valence/si2.amn")
    A1 = read_amn(artifact"Si2_valence/reference/binary/si2.amn")
    @test A ≈ A1

    tmpfile = tempname(; cleanup=true)
    write_amn(tmpfile, A; binary=true)
    A2 = read_amn(tmpfile)
    @test A ≈ A2
end
