@testset "read/write w90 band dat" begin
    band = read_w90_band(artifact"Si2/si2")

    outdir = mktempdir(; cleanup=true)
    outprefix = joinpath(outdir, "si2")

    write_w90_band(outprefix; band...)
    band2 = read_w90_band(outprefix)
    @test band == band2
end
