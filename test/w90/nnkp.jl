@testitem "read nnkp" begin
    using LazyArtifacts
    nnkp = read_nnkp(artifact"Si2_valence/outputs/Si2_valence.nnkp")

    WRITE_TOML = false
    WRITE_TOML && write_nnkp("/tmp/Si2_valence.nnkp.toml", nnkp, WannierIO.W90InputToml())

    test_data = read_nnkp(artifact"Si2_valence/outputs/Si2_valence.nnkp.toml")
    # make their keys unordered for comparison
    @test Dict(pairs(nnkp)) == Dict(pairs(test_data))
end

@testitem "read/write nnkp" begin
    using LazyArtifacts
    nnkp = read_nnkp(artifact"Si2_valence/outputs/Si2_valence.nnkp")
    tmpfile = tempname(; cleanup = true)
    write_nnkp(tmpfile, nnkp)

    nnkp2 = read_nnkp(tmpfile)
    @test nnkp == nnkp2
end

@testitem "read/write nnkp toml" begin
    using LazyArtifacts
    # Note that this requires https://github.com/JuliaLang/julia/pull/57584
    if VERSION > v"1.11.4"
        nnkp = read_nnkp(artifact"Si2_valence/outputs/Si2_valence.nnkp.toml")

        tmpfile = tempname(; cleanup = true)
        write_nnkp(tmpfile, nnkp, WannierIO.W90InputToml())

        nnkp2 = read_nnkp(tmpfile)
        @test nnkp == nnkp2
    end
end

@testitem "read spinor_projections" begin
    using LazyArtifacts
    nnkp = read_nnkp(artifact"Fe_soc/outputs/Fe.nnkp")

    @test haskey(nnkp, "spinor_projections")
    @test !haskey(nnkp, "projections")
    sp = nnkp["spinor_projections"]
    @test sp isa AbstractVector{<:WannierIO.SpinorHydrogenOrbital}
    @test length(sp) == 16
    # spin alternates +1 / -1 for consecutive up/down partners
    @test sp[1].spin == 1
    @test sp[2].spin == -1
    @test sp[1].spin_qaxis == [0.0, 0.0, 1.0]
end

@testitem "read/write spinor_projections" begin
    using LazyArtifacts
    nnkp = read_nnkp(artifact"Fe_soc/outputs/Fe.nnkp")
    tmpfile = tempname(; cleanup = true)
    write_nnkp(tmpfile, nnkp)

    nnkp2 = read_nnkp(tmpfile)
    @test nnkp == nnkp2
end

@testitem "read/write spinor_projections toml" begin
    using LazyArtifacts
    # Note that this requires https://github.com/JuliaLang/julia/pull/57584
    if VERSION > v"1.11.4"
        nnkp = read_nnkp(artifact"Fe_soc/outputs/Fe.nnkp")

        tmpfile = tempname(; cleanup = true)
        write_nnkp(tmpfile, nnkp, WannierIO.W90InputToml())

        nnkp2 = read_nnkp(tmpfile)
        @test nnkp == nnkp2
    end
end

@testitem "read auto_projections" begin
    using LazyArtifacts
    nnkp = read_nnkp(artifact"SnSe2/outputs/SnSe2.nnkp")
    @test nnkp["auto_projections"] == 12

    tmpfile = tempname(; cleanup = true)
    write_nnkp(tmpfile, nnkp, WannierIO.W90InputToml())
    nnkp2 = read_nnkp(tmpfile)
    @test nnkp["auto_projections"] == nnkp2["auto_projections"]
end
