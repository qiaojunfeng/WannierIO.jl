
@testset "read nnkp" begin
    toml_path = joinpath(@__DIR__, "test_data/nnkp.toml")

    nnkp = read_nnkp(joinpath(FIXTURE_PATH, "si2.nnkp"))

    WRITE_TOML = false
    WRITE_TOML && WannierIO._write_nnkp_toml(toml_path; nnkp...)

    test_data = WannierIO._read_nnkp_toml(toml_path)

    # make their keys unordered for comparison
    @test pairs(nnkp) == pairs(test_data)
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

@testset "read/write nnkp toml" begin
    toml_path = joinpath(@__DIR__, "test_data/nnkp.toml")
    nnkp = WannierIO._read_nnkp_toml(toml_path)

    tmpfile = tempname(; cleanup=true)
    WannierIO._write_nnkp_toml(tmpfile; nnkp...)

    nnkp2 = WannierIO._read_nnkp_toml(tmpfile)

    @test nnkp == nnkp2
end
