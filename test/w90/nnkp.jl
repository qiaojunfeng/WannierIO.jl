@testitem "read nnkp" begin
    using LazyArtifacts
    nnkp = read_nnkp(artifact"Si2_valence/reference/Si2_valence.nnkp")

    WRITE_TOML = false
    WRITE_TOML && write_nnkp("/tmp/Si2_valence.nnkp.toml"; toml=true, nnkp...)

    test_data = read_nnkp(artifact"Si2_valence/reference/Si2_valence.nnkp.toml")
    # make their keys unordered for comparison
    @test pairs(nnkp) == pairs(test_data)
end

@testitem "read/write nnkp" begin
    using LazyArtifacts
    nnkp = read_nnkp(artifact"Si2_valence/reference/Si2_valence.nnkp")
    tmpfile = tempname(; cleanup=true)
    n_wann = 4
    write_nnkp(tmpfile; nnkp..., n_wann)

    nnkp2 = read_nnkp(tmpfile)
    @test pairs(nnkp) == pairs(nnkp2)
end

@testitem "read/write nnkp toml" begin
    using LazyArtifacts
    nnkp = read_nnkp(artifact"Si2_valence/reference/Si2_valence.nnkp.toml")

    tmpfile = tempname(; cleanup=true)
    write_nnkp(tmpfile; toml=true, nnkp...)

    nnkp2 = read_nnkp(tmpfile)
    @test pairs(nnkp) == pairs(nnkp2)
end

@testitem "read auto_projections" begin
    using LazyArtifacts
    nnkp = read_nnkp(artifact"SnSe2/reference/SnSe2.nnkp")
    @test nnkp.auto_projections == 12

    tmpfile = tempname(; cleanup=true)
    write_nnkp(tmpfile; toml=true, nnkp...)
    nnkp2 = read_nnkp(tmpfile)
    @test nnkp.auto_projections == nnkp2.auto_projections
end
