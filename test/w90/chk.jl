@testset "read chk" begin
    chk = read_chk(artifact"Si2_valence/reference/si2.chk.fmt")

    @test chk.n_wann == 4
    @test chk.n_bands == 4
    @test chk.n_kpts == 8
    @test chk.n_bvecs == 8

    # make sure we read the lattice as column-major
    win = read_win(artifact"Si2_valence/si2.win")
    @test chk.lattice ≈ win.unit_cell_cart
end

@testset "read/write chk" begin
    chk = read_chk(artifact"Si2_valence/reference/si2.chk.fmt")

    tmpfile = tempname(; cleanup=true)
    write_chk(tmpfile, chk)
    chk2 = read_chk(tmpfile)
    @test chk ≈ chk2
end

@testset "read/write chk binary" begin
    chk = read_chk(artifact"Si2_valence/reference/si2.chk.fmt")
    chk1 = read_chk(artifact"Si2_valence/reference/binary/si2.chk")
    @test chk ≈ chk1

    tmpfile = tempname(; cleanup=true)
    write_chk(tmpfile, chk; binary=true)
    chk2 = read_chk(tmpfile)
    @test chk ≈ chk2
end

@testset "read chk disentanglement" begin
    chk = read_chk(artifact"Si2/reference/si2.chk.fmt")

    @test chk.n_wann == 4
    @test chk.n_bands == 12
    @test chk.n_kpts == 8
    @test chk.n_bvecs == 8

    # make sure we read the lattice as column-major
    win = read_win(artifact"Si2/si2.win")
    @test chk.lattice ≈ win.unit_cell_cart
end

@testset "read/write chk disentanglement" begin
    chk = read_chk(artifact"Si2/reference/si2.chk.fmt")

    tmpfile = tempname(; cleanup=true)
    write_chk(tmpfile, chk)
    chk2 = read_chk(tmpfile)
    @test chk ≈ chk2
end

@testset "read/write chk disentanglement binary" begin
    chk = read_chk(artifact"Si2/reference/si2.chk.fmt")
    chk1 = read_chk(artifact"Si2/reference/binary/si2.chk")
    @test chk ≈ chk1

    tmpfile = tempname(; cleanup=true)
    write_chk(tmpfile, chk; binary=true)
    chk2 = read_chk(tmpfile)
    @test chk ≈ chk2
end

@testset "get_Udis" begin
    chk = read_chk(artifact"Si2/reference/si2.chk")
    Udis = get_Udis(chk)
    Udis_ref = read_amn(artifact"Si2/reference/si2.chk_Udis.amn")
    @test Udis ≈ Udis_ref
end

@testset "get_U" begin
    chk = read_chk(artifact"Si2/reference/si2.chk")
    U = get_U(chk)
    U_ref = read_amn(artifact"Si2/reference/si2.chk_U.amn")
    @test U ≈ U_ref
end
