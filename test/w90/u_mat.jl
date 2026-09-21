@testitem "read/write u.mat" begin
    using LazyArtifacts
    chk = read_chk(artifact"Si2/outputs/Si2.chk.fmt")

    Udismat = WannierIO.read_u_mat(artifact"Si2/outputs/Si2_u.mat")
    @test Udismat.U ≈ chk.Uml
    @test Udismat.kpoints ≈ chk.kpoints

    tmpfile = tempname(; cleanup = true)
    WannierIO.write_u_mat(tmpfile, chk.Uml, chk.kpoints)
    Udismat2 = WannierIO.read_u_mat(tmpfile)
    @test Udismat2.U ≈ chk.Uml
    @test Udismat2.kpoints ≈ chk.kpoints
    # wannier90 precision loses ~1e-10; 16 digits round-trips a Float64 gauge
    @test maximum(abs, Udismat2.U - chk.Uml) > 1.0e-12
    WannierIO.write_u_mat(tmpfile, chk.Uml, chk.kpoints; digits = 16)
    Udismat3 = WannierIO.read_u_mat(tmpfile)
    @test maximum(abs, Udismat3.U - chk.Uml) < 1.0e-14
    @test Udismat3.kpoints ≈ chk.kpoints atol = 1.0e-15
    @test_throws ArgumentError WannierIO.write_u_mat(tmpfile, chk.Uml, chk.kpoints; digits = 0)
end

@testitem "read/write u_dis.mat" begin
    using LazyArtifacts
    chk = read_chk(artifact"Si2/outputs/Si2.chk.fmt")

    Udismat = WannierIO.read_u_mat(artifact"Si2/outputs/Si2_u_dis.mat")
    # do not use `gauge_matrices_dis` since it sorts the band indices, here we want to
    # compare the raw data
    @test Udismat.U ≈ chk.Udis
    @test Udismat.kpoints ≈ chk.kpoints

    tmpfile = tempname(; cleanup = true)
    WannierIO.write_u_mat(tmpfile, chk.Udis, chk.kpoints)
    Udismat2 = WannierIO.read_u_mat(tmpfile)
    @test Udismat2.U ≈ chk.Udis
    @test Udismat2.kpoints ≈ chk.kpoints
end
