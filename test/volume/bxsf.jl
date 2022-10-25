
@testset "read/write bxsf" begin
    bxsf = read_bxsf(joinpath(FIXTURE_PATH, "volume/copper.bxsf"))

    @test bxsf.fermi_energy ≈ 12.2102
    @test bxsf.origin ≈ [0.0, 0.0, 0.0]
    @test bxsf.span_vectors ≈
        [
        -1.7404719250572223 -1.7404719250572223 1.7404719250572223
        1.7404719250572223 1.7404719250572223 1.7404719250572223
        -1.7404719250572223 1.7404719250572223 -1.7404719250572223
    ]'
    @test bxsf.X ≈ range(0, 1, 3)
    @test bxsf.Y ≈ bxsf.X
    @test bxsf.Z ≈ bxsf.X
    @test size(bxsf.E) == (7, 3, 3, 3)

    tmpfile = tempname(; cleanup=true)
    write_bxsf(tmpfile, bxsf.fermi_energy, bxsf.origin, bxsf.span_vectors, bxsf.E)
    bxsf2 = read_bxsf(tmpfile)

    @test bxsf.fermi_energy ≈ bxsf2.fermi_energy
    @test bxsf.origin ≈ bxsf2.origin
    @test bxsf.span_vectors ≈ bxsf2.span_vectors
    @test bxsf.X ≈ bxsf2.X
    @test bxsf.Y ≈ bxsf2.Y
    @test bxsf.Z ≈ bxsf2.Z
    @test bxsf.E ≈ bxsf2.E
end
