
@testset "read nnkp" begin
    test_data = YAML.load_file(String(@__DIR__) * "/test_data/nnkp.yaml")

    nnkp = read_nnkp(joinpath(FIXTURE_PATH, "si2.nnkp"))

    # Convert type so YAML can write it.
    kpb_b = nnkp.kpb_b
    dict = Dict(
        "recip_lattice" => [[nnkp.recip_lattice[i,j] for i = 1:3] for j = 1:3],
        "kpoints" => nnkp.kpoints,
        "kpb_k" => nnkp.kpb_k,
        "kpb_b" => nnkp.kpb_b
    )


    # YAML.write_file(String(@__DIR__) * "/test_data/nnkp.yaml", dict)

    for (key, value) in dict
        @test value ≈ test_data[key]
    end
end

@testset "read/write nnkp" begin
    nnkp = read_nnkp(joinpath(FIXTURE_PATH, "si2.nnkp"))
    tmpfile = tempname(; cleanup=true)
    n_wann = 8
    write_nnkp(tmpfile, nnkp.recip_lattice, nnkp.kpoints, nnkp.kpb_k, nnkp.kpb_b, n_wann)

    nnkp2 = read_nnkp(tmpfile)
    @test nnkp.recip_lattice ≈ nnkp2.recip_lattice
    @test nnkp.kpoints ≈ nnkp2.kpoints
    @test nnkp.kpb_k ≈ nnkp2.kpb_k
    @test nnkp.kpb_b ≈ nnkp2.kpb_b
end
