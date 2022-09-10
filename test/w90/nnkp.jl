
@testset "read nnkp" begin
    test_data = YAML.load_file(String(@__DIR__) * "/test_data/nnkp.yaml")

    bvectors = read_nnkp(joinpath(FIXTURE_PATH, "silicon/silicon.nnkp"))

    # Convert type so YAML can write it.
    kpb_b = bvectors.kpb_b
    dict = Dict(
        "recip_lattice" => mat2vec(bvectors.recip_lattice),
        "kpoints" => mat2vec(bvectors.kpoints),
        "bvectors" => mat2vec(bvectors.bvectors),
        "kpb_k" => mat2vec(bvectors.kpb_k),
        "kpb_b" => [mat2vec(kpb_b[:, :, ik]) for ik in axes(kpb_b, 3)],
    )

    # YAML.write_file(String(@__DIR__) * "/test_data/nnkp.yaml", dict)

    for (key, value) in dict
        @test value ≈ test_data[key]
    end
end

@testset "read/write nnkp" begin
    bvectors = read_nnkp(joinpath(FIXTURE_PATH, "silicon/silicon.nnkp"))
    tmpfile = tempname(; cleanup=true)
    n_wann = 8
    write_nnkp(tmpfile, bvectors, n_wann)

    bvectors2 = read_nnkp(tmpfile)
    @test bvectors.recip_lattice ≈ bvectors2.recip_lattice
    @test bvectors.kpoints ≈ bvectors2.kpoints
    @test bvectors.bvectors ≈ bvectors2.bvectors
    @test bvectors.kpb_k ≈ bvectors2.kpb_k
    @test bvectors.kpb_b ≈ bvectors2.kpb_b
end
