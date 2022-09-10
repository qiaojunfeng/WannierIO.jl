
@testset "read wsvec" begin
    Rvecs = Wannier.read_w90_wsvec(
        joinpath(FIXTURE_PATH, "valence/band/ws/silicon_wsvec.dat")
    )
    # just some simple tests
    @test Rvecs.R[:, 1] == [-3, 1, 1]

    Rvecs = Wannier.read_w90_wsvec(
        joinpath(FIXTURE_PATH, "valence/band/mdrs/silicon_wsvec.dat")
    )
    @test Rvecs.R[:, 1] == [-3, 1, 1]
end

@testset "read tb" begin
    Rvecs, H, positions = Wannier.read_w90_tb(
        joinpath(FIXTURE_PATH, "valence/band/ws/silicon")
    )
    # just some simple tests
    R1 = [-3, 1, 1]
    @test Rvecs.R[:, 1] == R1
    H111 = 0.51893360E-02 + im * -0.29716277E-02
    @test H[1, 1, 1] ≈ H111
    P111end = 0.24832468E-03 + im * -0.21054981E-03
    @test positions[1, 1, 1, end] ≈ P111end

    Rvecs, H, positions = Wannier.read_w90_tb(
        joinpath(FIXTURE_PATH, "valence/band/mdrs/silicon")
    )
    @test Rvecs.R[:, 1] == R1
    @test Rvecs.T[1, 1, 1] == [0 4 4 4; 0 -4 0 0; 0 0 -4 0]
    @test Rvecs.Nᵀ[1, 1, 1] == 4
    @test H[1, 1, 1] ≈ H111
    @test positions[1, 1, 1, end] ≈ P111end
end
