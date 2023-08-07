@testitem "read/write spn" begin
    using LazyArtifacts
    Sx, Sy, Sz = read_spn(artifact"Fe_soc/test/Fe.spn.fmt")

    ref_Sx_end = ComplexF64[
        -1.045982416783206e-9+0.0im 1.8309573615290023e-9+1.8135039310434643e-9im -0.9143596596228599+0.39269798727353im 2.2580234889560057e-5+0.00012322780868176796im
        1.8309573615290023e-9-1.8135039310434643e-9im 5.005565130696888e-9+0.0im -6.863936845424341e-5+0.00010788172874406267im 0.6996261651438792+0.7076597058385461im
        -0.9143596596228599-0.39269798727353im -6.863936845424341e-5-0.00010788172874406267im 8.544459741957264e-9+0.0im -1.1837183322351442e-8+7.277568287854686e-10im
        2.2580234889560057e-5-0.00012322780868176796im 0.6996261651438792-0.7076597058385461im -1.1837183322351442e-8-7.277568287854686e-10im -1.3997437777976164e-8+0.0im
    ]

    @test length(Sx) == length(Sy) == length(Sz) == 8
    @test size(Sx[1]) == size(Sy[1]) == size(Sz[1]) == (4, 4)
    @test Sx[end] ≈ ref_Sx_end

    tmpfile = tempname(; cleanup=true)
    write_spn(tmpfile, Sx, Sy, Sz; binary=false)

    Sx1, Sy1, Sz1 = read_spn(tmpfile)
    @test Sx ≈ Sx1
    @test Sy ≈ Sy1
    @test Sz ≈ Sz1
end

@testitem "read/write spn binary" begin
    using LazyArtifacts
    Sx, Sy, Sz = read_spn(artifact"Fe_soc/test/Fe.spn")
    Sx1, Sy1, Sz1 = read_spn(artifact"Fe_soc/test/Fe.spn.fmt")
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
