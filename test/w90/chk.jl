
@testset "read chk" begin
    chk = read_chk(joinpath(FIXTURE_PATH, "silicon/silicon.chk.fmt"))

    @test chk.n_wann == 8
    @test chk.n_bands == 12
end

@testset "read/write chk" begin
    chk = read_chk(joinpath(FIXTURE_PATH, "silicon/silicon.chk.fmt"))

    win = read_win(joinpath(FIXTURE_PATH, "silicon/silicon.win"))
    @test chk.lattice ≈ win.unit_cell

    tmpfile = tempname(; cleanup=true)
    write_chk(tmpfile, chk)
    chk2 = read_chk(tmpfile)

    @test chk.header == chk2.header
    @test chk.n_bands == chk2.n_bands
    @test chk.n_exclude_bands == chk2.n_exclude_bands
    @test chk.exclude_bands == chk2.exclude_bands
    @test chk.lattice ≈ chk2.lattice
    @test chk.recip_lattice ≈ chk2.recip_lattice
    @test chk.n_kpts ≈ chk2.n_kpts
    @test chk.kgrid == chk2.kgrid
    @test chk.kpoints ≈ chk2.kpoints
    @test chk.n_bvecs == chk2.n_bvecs
    @test chk.n_wann == chk2.n_wann
    @test chk.checkpoint == chk2.checkpoint
    @test chk.have_disentangled == chk2.have_disentangled
    @test chk.ΩI ≈ chk2.ΩI
    @test chk.dis_bands == chk2.dis_bands
    @test chk.Uᵈ ≈ chk2.Uᵈ
    @test chk.U ≈ chk2.U
    @test chk.M ≈ chk2.M
    @test chk.r ≈ chk2.r
    @test chk.ω ≈ chk2.ω
    @test chk.n_dis == chk2.n_dis
end
