@testset "read/write w90 band dat" begin
    band = read_w90_band(artifact"Si2_valence/reference/mdrs/si2")

    @test length(band.kpoints) == 511
    @test length(band.kpoints[1]) == 3
    @test band.kpoints[511] == [0.5, 0.0, 0.5]
    @test band.symm_idx == [1, 101, 136, 197, 303, 390, 461, 511]
    @test band.symm_label == ["G", "X", "U", "K", "G", "L", "W", "X"]
    @test length(band.x) == 511
    @test band.x[511] == 5.9004323

    outdir = mktempdir(; cleanup=true)
    outprefix = joinpath(outdir, "si2")

    write_w90_band(outprefix; band...)
    band2 = read_w90_band(outprefix)
    @test band == band2
end
