
@testset "read/write mmn" begin
    M, kpb_k, kpb_b = read_mmn(joinpath(FIXTURE_PATH, "formatted/si2.mmn"))

    tmpfile = tempname(; cleanup=true)
    write_mmn(tmpfile, M, kpb_k, kpb_b)
    M2, kpb_k2, kpb_b2 = read_mmn(tmpfile)

    @test M ≈ M2
    @test kpb_k ≈ kpb_k2
    @test kpb_b ≈ kpb_b2
end

@testset "read/write mmn binary" begin
    M, kpb_k, kpb_b = read_mmn(joinpath(FIXTURE_PATH, "formatted/si2.mmn"))
    M1, kpb_k1, kpb_b1 = read_mmn(joinpath(FIXTURE_PATH, "unformatted/si2.mmn"))
    @test M ≈ M1
    @test kpb_k ≈ kpb_k1
    @test kpb_b ≈ kpb_b1

    tmpfile = tempname(; cleanup=true)
    write_mmn(tmpfile, M, kpb_k, kpb_b; binary=true)
    M2, kpb_k2, kpb_b2 = read_mmn(tmpfile)

    @test M ≈ M2
    @test kpb_k ≈ kpb_k2
    @test kpb_b ≈ kpb_b2
end
