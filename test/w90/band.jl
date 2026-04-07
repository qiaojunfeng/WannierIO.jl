@testitem "read/write w90 band dat" begin
    using LazyArtifacts
    band = read_w90_band(artifact"Si2_valence/outputs/MDRS/Si2_valence")

    @test length(band.kpoints) == 511
    @test length(band.kpoints[1]) == 3
    @test band.kpoints[511] == [0.5, 0.0, 0.5]
    @test band.symm_point_indices == [1, 101, 136, 197, 303, 390, 461, 511]
    @test band.symm_point_labels == ["G", "X", "U", "K", "G", "L", "W", "X"]
    @test length(band.x) == 511
    @test band.x[511] == 5.9004323

    outdir = mktempdir(; cleanup = true)
    outprefix = joinpath(outdir, "Si2_valence")

    write_w90_band(outprefix; band...)
    band2 = read_w90_band(outprefix)
    @test band == band2
end

@testitem "read/write w90 band" begin
    using LazyArtifacts

    win = read_win(artifact"Si2_valence/Si2_valence.win")
    recip_lattice = WannierIO.reciprocal_lattice(win["unit_cell_cart"])
    kpath, eigenvalues = read_w90_band(
        artifact"Si2_valence/outputs/MDRS/Si2_valence", recip_lattice
    )

    outdir = mktempdir(; cleanup = true)
    outprefix = joinpath(outdir, "Si2_valence")
    write_w90_band(outprefix, kpath, eigenvalues)
    kpath2, eigenvalues2 = read_w90_band(outprefix, recip_lattice)

    @test kpath ≈ kpath2
    @test eigenvalues ≈ eigenvalues2
end

@testitem "read/write_w90_band_kpt_labelinfo" begin
    using LazyArtifacts

    win = read_win(artifact"Si2_valence/Si2_valence.win")
    recip_lattice = WannierIO.reciprocal_lattice(win["unit_cell_cart"])
    kpath, eigenvalues = read_w90_band(
        artifact"Si2_valence/outputs/MDRS/Si2_valence", recip_lattice
    )

    outdir = mktempdir(; cleanup = true)
    outprefix = joinpath(outdir, "Si2_valence")
    WannierIO.write_w90_band_kpt_labelinfo(outprefix, kpath)
    kpath2 = WannierIO.read_w90_band_kpt_labelinfo(outprefix, recip_lattice)

    @test kpath ≈ kpath2
end
