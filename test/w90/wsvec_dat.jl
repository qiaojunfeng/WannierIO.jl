@testitem "read wsvec WS" begin
    using LazyArtifacts
    wsvec = read_w90_wsvec_dat(artifact"Si2_valence/outputs/WS/Si2_valence_wsvec.dat")

    @assert wsvec.mdrs == false
    @test length(wsvec.Rvectors) == 279
    @test wsvec.Rvectors[1] == [-4, 0, 2]
    @test wsvec.n_wann == 4
end

@testitem "read wsvec MDRS" begin
    using WannierIO: Vec3
    using LazyArtifacts
    wsvec = read_w90_wsvec_dat(artifact"Si2_valence/outputs/MDRS/Si2_valence_wsvec.dat")

    @assert wsvec.mdrs == true
    @test length(wsvec.Rvectors) == 279
    @test wsvec.Rvectors[1] == [-4, 0, 2]
    @test size(wsvec.Tvectors) == (4, 4, 279)
    @test wsvec.Tvectors[1, 1, 1] == Vector{Vec3}([[0, 0, 0], [6, 0, -6], [6, 0, 0]])
    @test size(wsvec.Tdegens) == (4, 4, 279)
    @test wsvec.Tdegens[:, :, 1] == [3 1 1 1; 1 3 1 1; 2 2 3 2; 1 1 1 3]
    @test wsvec.n_wann == 4
end

@testitem "WsvecDat constructors" begin
    using WannierIO: Vec3

    Rvectors = [Vec3(0, 0, 0)]
    wsvec = WannierIO.WsvecDat("constructor test", Rvectors, 1)

    @test wsvec.header == "constructor test"
    @test wsvec.mdrs == false
    @test wsvec.Rvectors == Rvectors
    @test wsvec.n_wann == 1

    Tvectors = Array{Vector{Vec3{Int}}, 3}(undef, 1, 1, 1)
    Tvectors[1, 1, 1] = [Vec3(0, 0, 0)]
    Tdegens = zeros(Int, 1, 1, 1)
    Tdegens[1, 1, 1] = 1
    wsvec_mdrs = WannierIO.WsvecDat("constructor test", Rvectors, Tvectors, Tdegens)

    @test wsvec_mdrs.header == "constructor test"
    @test wsvec_mdrs.mdrs == true
    @test wsvec_mdrs.Rvectors == Rvectors
    @test wsvec_mdrs.Tvectors == Tvectors
    @test wsvec_mdrs.Tdegens == Tdegens
    @test wsvec_mdrs.n_wann == 1
end

@testitem "write wsvec WS" begin
    using LazyArtifacts
    wsvec = read_w90_wsvec_dat(artifact"Si2_valence/outputs/WS/Si2_valence_wsvec.dat")

    tmpfile = tempname(; cleanup = true)
    write_w90_wsvec_dat(tmpfile, wsvec)
    wsvec2 = read_w90_wsvec_dat(tmpfile)

    @test propertynames(wsvec) == propertynames(wsvec2)
    for name in propertynames(wsvec)
        name == :header && continue
        @test getfield(wsvec2, name) == getfield(wsvec, name)
    end
end

@testitem "write wsvec MDRS" begin
    using WannierIO: Vec3
    using LazyArtifacts
    wsvec = read_w90_wsvec_dat(artifact"Si2_valence/outputs/MDRS/Si2_valence_wsvec.dat")

    tmpfile = tempname(; cleanup = true)
    write_w90_wsvec_dat(tmpfile, wsvec)
    wsvec2 = read_w90_wsvec_dat(tmpfile)

    @test propertynames(wsvec) == propertynames(wsvec2)
    for name in propertynames(wsvec)
        name == :header && continue
        @test getfield(wsvec2, name) == getfield(wsvec, name)
    end

    # n_wann is optional
    write_w90_wsvec_dat(tmpfile, wsvec)
    wsvec2 = read_w90_wsvec_dat(tmpfile)

    @test propertynames(wsvec) == propertynames(wsvec2)
    for name in propertynames(wsvec)
        name == :header && continue
        @test getfield(wsvec2, name) == getfield(wsvec, name)
    end
end
