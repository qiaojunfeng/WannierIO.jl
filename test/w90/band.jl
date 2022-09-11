@testset "read/write w90 band dat" begin
    band = read_w90_band(joinpath(FIXTURE_PATH, "si2"))

    outdir = mktempdir(; cleanup=true)
    outseedname = joinpath(outdir, "si2")

    write_w90_band(
        outseedname, band.kpoints, band.E, band.x, band.symm_idx, band.symm_label
    )

    band2 = read_w90_band(outseedname)

    @test band.kpoints ≈ band2.kpoints
    @test band.E ≈ band2.E
    @test band.x ≈ band2.x
    @test band.symm_idx == band2.symm_idx
    @test band.symm_label == band2.symm_label
end
