@testitem "read/write w90 band dat" begin
    using LazyArtifacts
    band = read_w90_band(artifact"Si2_valence/outputs/MDRS/Si2_valence")

    @test length(band.kpoints) == 511
    @test length(band.kpoints[1]) == 3
    @test band.kpoints[511] == [0.5, 0.0, 0.5]
    @test band.symm_point_indices == [1, 101, 136, 197, 303, 390, 461, 511]
    @test band.symm_point_labels == ["G", "X", "U", "K", "G", "L", "W", "X"]
    @test length(band.x) == 511
    @test band.x[511] == 5.9004323

    outdir = mktempdir(; cleanup=true)
    outprefix = joinpath(outdir, "Si2_valence")

    write_w90_band(outprefix; band...)
    band2 = read_w90_band(outprefix)
    @test band == band2
end
