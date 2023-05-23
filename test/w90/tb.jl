using WannierIO: Vec3
@testset "read wsvec" begin
    mdrs, wsvec = read_w90_wsvec(joinpath(FIXTURE_PATH, "ws/si2_wsvec.dat"))
    # just some simple tests
    @assert mdrs == false
    @test wsvec.R[1] == [-1, -1, 1]

    mdrs, wsvec = read_w90_wsvec(joinpath(FIXTURE_PATH, "mdrs/si2_wsvec.dat"))
    @assert mdrs == true
    @test wsvec.R[1] == [-1, -1, 1]
    @test wsvec.T[1][1, 1] == [Vec3([0 0 0 2 2 2; 0 2 2 0 0 2; 0 -2 0 -2 0 -2][:, i]) for i = 1:6] 
    @test wsvec.Nᵀ[1][1, 1] == 6
end

@testset "read tb" begin
    tbdat = read_w90_tbdat(joinpath(FIXTURE_PATH, "ws/si2_tb.dat"))
    # just some simple tests
    R1 = [-1, -1, 1]
    @test tbdat.R[1] == R1
    H111 = 0.12534576E-02 + im * -0.75115190E-17
    @test tbdat.H[1][1, 1] ≈ H111
    r111end = 0.17347235E-17 + im * -0.92553612E-23
    @test tbdat.Rx[end][1, 1] ≈ r111end

    tbdat = read_w90_tbdat(joinpath(FIXTURE_PATH, "mdrs/si2_tb.dat"))
    @test tbdat.R[1] == R1
    @test tbdat.H[1][1, 1] ≈ H111
    @test tbdat.Rx[end][1, 1] ≈ r111end
end
