@testitem "read wsvec WS" begin
    using LazyArtifacts
    wsvec = read_w90_wsvec(artifact"Si2_valence/outputs/WS/Si2_valence_wsvec.dat")

    @assert wsvec.mdrs == false
    @test length(wsvec.Rvectors) == 279
    @test wsvec.Rvectors[1] == [-4, 0, 2]
    @test wsvec.n_wann == 4
end

@testitem "read wsvec MDRS" begin
    using WannierIO: Vec3
    using LazyArtifacts
    wsvec = read_w90_wsvec(artifact"Si2_valence/outputs/MDRS/Si2_valence_wsvec.dat")

    @assert wsvec.mdrs == true
    @test length(wsvec.Rvectors) == 279
    @test wsvec.Rvectors[1] == [-4, 0, 2]
    @test length(wsvec.Tvectors) == 279
    @test wsvec.Tvectors[1][1, 1] == Vector{Vec3}([[0, 0, 0], [6, 0, -6], [6, 0, 0]])
    @test length(wsvec.Tdegens) == 279
    @test wsvec.Tdegens[1] == [3 1 1 1; 1 3 1 1; 2 2 3 2; 1 1 1 3]
    @test wsvec.n_wann == 4
end

@testitem "write wsvec WS" begin
    using LazyArtifacts
    wsvec = read_w90_wsvec(artifact"Si2_valence/outputs/WS/Si2_valence_wsvec.dat")

    tmpfile = tempname(; cleanup=true)
    write_w90_wsvec(tmpfile; wsvec.Rvectors, wsvec.n_wann)
    wsvec2 = read_w90_wsvec(tmpfile)

    @test keys(wsvec) == keys(wsvec2)
    for (k, v) in pairs(wsvec)
        k == :header && continue
        @test wsvec2[k] == v
    end
end

@testitem "write wsvec MDRS" begin
    using WannierIO: Vec3
    using LazyArtifacts
    wsvec = read_w90_wsvec(artifact"Si2_valence/outputs/MDRS/Si2_valence_wsvec.dat")

    tmpfile = tempname(; cleanup=true)
    write_w90_wsvec(tmpfile; wsvec.Rvectors, wsvec.Tvectors, wsvec.Tdegens)
    wsvec2 = read_w90_wsvec(tmpfile)

    @test keys(wsvec) == keys(wsvec2)
    for (k, v) in pairs(wsvec)
        k == :header && continue
        @test wsvec2[k] == v
    end

    # n_wann is optional
    write_w90_wsvec(tmpfile; wsvec.Rvectors, wsvec.Tvectors, wsvec.Tdegens, wsvec.n_wann)
    wsvec2 = read_w90_wsvec(tmpfile)

    @test keys(wsvec) == keys(wsvec2)
    for (k, v) in pairs(wsvec)
        k == :header && continue
        @test wsvec2[k] == v
    end
end
