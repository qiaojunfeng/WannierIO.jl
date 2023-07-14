
@testset "read nnkp" begin
    toml_path = joinpath(@__DIR__, "test_data/nnkp.toml")
    nnkp = read_nnkp(joinpath(FIXTURE_PATH, "si2.nnkp"))

    WRITE_TOML = false
    WRITE_TOML && write_nnkp(toml_path; toml=true, nnkp...)

    test_data = read_nnkp(toml_path)
    # make their keys unordered for comparison
    @test pairs(nnkp) == pairs(test_data)
end

@testset "read/write nnkp" begin
    nnkp = read_nnkp(joinpath(FIXTURE_PATH, "si2.nnkp"))

    tmpfile = tempname(; cleanup=true)
    n_wann = 8
    write_nnkp(tmpfile; nnkp..., n_wann)

    nnkp2 = read_nnkp(tmpfile)
    @test nnkp == nnkp2
end

@testset "read/write nnkp toml" begin
    toml_path = joinpath(@__DIR__, "test_data/nnkp.toml")
    nnkp = read_nnkp(toml_path)

    tmpfile = tempname(; cleanup=true)
    write_nnkp(tmpfile; toml=true, nnkp...)

    nnkp2 = read_nnkp(tmpfile)
    @test nnkp == nnkp2
end
