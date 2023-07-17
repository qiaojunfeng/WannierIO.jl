@testset "read nnkp" begin
    nnkp = read_nnkp(artifact"Si2_valence/reference/Si2_valence.nnkp")

    WRITE_TOML = false
    WRITE_TOML && write_nnkp("/tmp/Si2_valence.nnkp.toml"; toml=true, nnkp...)

    test_data = read_nnkp(artifact"Si2_valence/reference/Si2_valence.nnkp.toml")
    # make their keys unordered for comparison
    @test pairs(nnkp) == pairs(test_data)
end

@testset "read/write nnkp" begin
    nnkp = read_nnkp(artifact"Si2_valence/reference/Si2_valence.nnkp")
    tmpfile = tempname(; cleanup=true)
    n_wann = 4
    write_nnkp(tmpfile; nnkp..., n_wann)

    nnkp2 = read_nnkp(tmpfile)
    @test pairs(nnkp) == pairs(nnkp2)
end

@testset "read/write nnkp toml" begin
    nnkp = read_nnkp(artifact"Si2_valence/reference/Si2_valence.nnkp.toml")

    tmpfile = tempname(; cleanup=true)
    write_nnkp(tmpfile; toml=true, nnkp...)

    nnkp2 = read_nnkp(tmpfile)
    @test pairs(nnkp) == pairs(nnkp2)
end
