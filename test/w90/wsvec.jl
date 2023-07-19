@testitem "read wsvec ws" begin
    using LazyArtifacts
    wsvec = read_w90_wsvec(artifact"Si2_valence/reference/ws/Si2_valence_wsvec.dat")

    @assert wsvec.mdrs == false
    @test length(wsvec.Rvectors) == 279
    @test wsvec.Rvectors[1] == [-4, 0, 2]
end

@testitem "read wsvec mdrs" begin
    using WannierIO: Vec3
    using LazyArtifacts
    wsvec = read_w90_wsvec(artifact"Si2_valence/reference/mdrs/Si2_valence_wsvec.dat")

    @assert wsvec.mdrs == true
    @test length(wsvec.Rvectors) == 279
    @test wsvec.Rvectors[1] == [-4, 0, 2]
    @test length(wsvec.Tvectors) == 279
    @test wsvec.Tvectors[1][1, 1] == Vector{Vec3}([[0, 0, 0], [6, 0, -6], [6, 0, 0]])
    @test length(wsvec.Tdegens) == 279
    @test wsvec.Tdegens[1] == [3 1 1 1; 1 3 1 1; 2 2 3 2; 1 1 1 3]
end
