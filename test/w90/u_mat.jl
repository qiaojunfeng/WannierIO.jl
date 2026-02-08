@testitem "read/write u.mat" begin
    using LazyArtifacts
    chk = read_chk(artifact"Si2/outputs/Si2.chk.fmt")

    Udismat = WannierIO.read_u_mat(artifact"Si2/outputs/Si2_u.mat")
    @test Udismat.U ≈ chk.Uml
    @test Udismat.kpoints ≈ chk.kpoints

    tmpfile = tempname(; cleanup=true)
    WannierIO.write_u_mat(tmpfile, chk.Uml, chk.kpoints)
    Udismat2 = WannierIO.read_u_mat(tmpfile)
    @test Udismat2.U ≈ chk.Uml
    @test Udismat2.kpoints ≈ chk.kpoints
end

@testitem "read/write u_dis.mat" begin
    using LazyArtifacts
    chk = read_chk(artifact"Si2/outputs/Si2.chk.fmt")

    Udismat = WannierIO.read_u_mat(artifact"Si2/outputs/Si2_u_dis.mat")
    # do not use `get_Udis` since it sorts the band indices, here we want to
    # compare the raw data
    @test Udismat.U ≈ chk.Udis
    @test Udismat.kpoints ≈ chk.kpoints

    tmpfile = tempname(; cleanup=true)
    WannierIO.write_u_mat(tmpfile, chk.Udis, chk.kpoints)
    Udismat2 = WannierIO.read_u_mat(tmpfile)
    @test Udismat2.U ≈ chk.Udis
    @test Udismat2.kpoints ≈ chk.kpoints
end
