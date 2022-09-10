
@testset "read/write mmn" begin
    M, kpb_k, kpb_b = read_mmn(joinpath(FIXTURE_PATH, "silicon/silicon.mmn"))

    tmpfile = tempname(; cleanup=true)

    write_mmn(tmpfile, M, kpb_k, kpb_b)

    M2, kpb_k2, kpb_b2 = read_mmn(tmpfile)

    @test M ≈ M2
    @test kpb_k ≈ kpb_k2
    @test kpb_b ≈ kpb_b2
end
