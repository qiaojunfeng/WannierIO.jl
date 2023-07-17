@testset "read/write bxsf" begin
    bxsf = read_bxsf(artifact"Cu/reference/Cu.bxsf")

    @test bxsf.fermi_energy ≈ 16.8985
    @test bxsf.origin ≈ [0.0, 0.0, 0.0]
    @test bxsf.span_vectors ≈
        [0.0 2.1769775 2.1769775; 2.1769775 0.0 2.1769775; 2.304768 2.304768 0.0]
    @test bxsf.X ≈ range(0, 1, 31)
    @test bxsf.Y ≈ bxsf.X
    @test bxsf.Z ≈ bxsf.X
    @test size(bxsf.E) == (7, 31, 31, 31)
    E1234 = 8.26105181
    @test bxsf.E[1, 2, 3, 4] ≈ E1234

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
