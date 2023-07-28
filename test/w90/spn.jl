@testitem "read/write spn" begin
    using LazyArtifacts
    Sx, Sy, Sz = read_spn(artifact"Fe/reference/Fe.spn.fmt")

    ref_Sx_end_3 = ComplexF64[
        1.227127937786662e-11-3.4992442140773507e-25im -0.028729507210608192+0.9982352627226101im -8.138264000334277e-10-3.415217248618749e-9im
        -0.028729507210608192-0.9982352627226101im 1.3574459202148256e-10-4.324777379222602e-26im -2.438410338009167e-11-4.515460432474117e-11im
        -8.138264000334277e-10+3.415217248618749e-9im -2.438410338009167e-11+4.515460432474117e-11im 4.9378004394971065e-12+2.220507462135635e-25im
    ]

    @test length(Sx) == length(Sy) == length(Sz) == 216
    @test size(Sx[1]) == size(Sy[1]) == size(Sz[1]) == (22, 22)
    @test Sx[end][1:3, 1:3] ≈ ref_Sx_end_3

    tmpfile = tempname(; cleanup=true)
    write_spn(tmpfile, Sx, Sy, Sz; binary=false)

    Sx1, Sy1, Sz1 = read_spn(tmpfile)
    @test Sx ≈ Sx1
    @test Sy ≈ Sy1
    @test Sz ≈ Sz1
end

@testitem "read/write spn binary" begin
    using LazyArtifacts
    Sx, Sy, Sz = read_spn(artifact"Fe/Fe.spn")
    Sx1, Sy1, Sz1 = read_spn(artifact"Fe/reference/Fe.spn.fmt")
    @test Sx ≈ Sx1
    @test Sy ≈ Sy1
    @test Sz ≈ Sz1

    tmpfile = tempname(; cleanup=true)
    write_spn(tmpfile, Sx, Sy, Sz; binary=true)
    Sx2, Sy2, Sz2 = read_spn(tmpfile)
    @test Sx ≈ Sx2
    @test Sy ≈ Sy2
    @test Sz ≈ Sz2
end
