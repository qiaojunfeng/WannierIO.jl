
@testset "read nnkp" begin
    nnkp = read_nnkp(artifact"Si2_valence/reference/si2.nnkp")

    WRITE_TOML = false
    WRITE_TOML && WannierIO._write_nnkp_toml("/tmp/si2.nnkp.toml"; nnkp...)

    test_data = WannierIO._read_nnkp_toml(artifact"Si2_valence/reference/si2.nnkp.toml")

    # make their keys unordered for comparison
    @test pairs(nnkp) == pairs(test_data)
end

@testset "read/write nnkp" begin
    nnkp = read_nnkp(artifact"Si2_valence/reference/si2.nnkp")
    tmpfile = tempname(; cleanup=true)
    n_wann = 4
    write_nnkp(tmpfile, nnkp.recip_lattice, nnkp.kpoints, nnkp.kpb_k, nnkp.kpb_b, n_wann)

    nnkp2 = read_nnkp(tmpfile)
    @test nnkp == nnkp2
end

@testset "read/write nnkp toml" begin
    nnkp = WannierIO._read_nnkp_toml(artifact"Si2_valence/reference/si2.nnkp.toml")

    tmpfile = tempname(; cleanup=true)
    write_nnkp(tmpfile; toml=true, nnkp...)

    nnkp2 = read_nnkp(tmpfile)
    @test nnkp == nnkp2
end
