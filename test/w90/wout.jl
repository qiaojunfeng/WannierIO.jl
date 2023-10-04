@testitem "read wout" begin
    using WannierIO: Vec3
    using LazyArtifacts
    wout = read_wout(artifact"Si2_valence/reference/Si2_valence.wout")

    ref_lattice = [
        0.000000 2.715265 2.715265
        2.715265 0.000000 2.715265
        2.715265 2.715265 0.000000
    ]
    ref_atom_labels = ["Si", "Si"]
    ref_atom_positions = Vec3[[0.00000, 0.00000, 0.00000], [0.25000, 0.25000, 0.25000]]
    ref_centers = Vec3[
        [0.678816, -0.678816, -0.678816],
        [-0.678816, -0.678816, 0.678816],
        [-0.678816, 0.678816, -0.678816],
        [0.678816, 0.678816, 0.678816],
    ]
    ref_spreads = [1.92917892, 1.92917900, 1.92917883, 1.92917887]

    @test wout.lattice ≈ ref_lattice
    @test wout.atom_labels == ref_atom_labels
    @test wout.atom_positions ≈ ref_atom_positions
    @test wout.centers ≈ ref_centers
    @test wout.spreads ≈ ref_spreads

    ref_ΩI = 7.153329184
    ref_ΩD = 0.000000000
    ref_ΩOD = 0.563386215
    ref_Ωtotal = 7.716715399
    @test wout.ΩI ≈ ref_ΩI
    @test wout.ΩD ≈ ref_ΩD
    @test wout.ΩOD ≈ ref_ΩOD
    @test wout.Ωtotal ≈ ref_Ωtotal
end

# test for non-symmetric lattice, with disentanglement
@testitem "read wout disentanglement" begin
    using LazyArtifacts
    wout = read_wout(artifact"Fe_soc/reference/Fe.wout")

    ref_lattice = [
        1.434996 -1.434996 -1.434996
        1.434996 1.434996 -1.434996
        1.434996 1.434996 1.434996
    ]
    ref_recip_lattice = [
        2.189269 -2.189269 0.0
        0.0 2.189269 -2.189269
        2.189269 0.0 2.189269
    ]
    ref_atom_labels = ["Fe"]
    ref_atom_positions = [[0.00000, 0.00000, 0.00000]]
    ref_centers = [
        [-0.001227, 0.761822, -1.210363],
        [-0.001132, 0.761153, 1.211403],
        [1.43403, 0.674118, -0.223106],
        [-1.436981, -0.000412, 0.95882],
        [-0.002088, -1.435412, 0.479468],
        [-0.001066, -0.761483, -1.210896],
        [0.00033, -0.000261, 0.000357],
        [8.5e-5, -0.019916, -0.01021],
        [-0.000258, -0.000766, 0.000475],
        [6.7e-5, -0.00072, -0.001011],
        [8.2e-5, 0.020905, 0.010305],
        [4.5e-5, 0.019902, -0.010005],
        [-0.000815, 0.00031, -0.000449],
        [0.000117, -0.020922, 0.010472],
        [-0.000259, 0.000986, 0.001254],
        [5.5e-5, 0.000493, -0.00047],
    ]
    ref_spreads = [
        1.11313978,
        1.11389149,
        1.11392684,
        1.14257541,
        1.14112348,
        1.11318311,
        0.40981481,
        0.44675647,
        0.48652472,
        0.44072505,
        0.45362382,
        0.44649383,
        0.4400766,
        0.45389206,
        0.48648695,
        0.44062982,
    ]

    @test wout.lattice ≈ ref_lattice
    @test wout.recip_lattice ≈ ref_recip_lattice
    @test wout.atom_labels == ref_atom_labels
    @test wout.atom_positions ≈ ref_atom_positions
    @test wout.centers ≈ ref_centers
    @test wout.spreads ≈ ref_spreads

    ref_ΩI = 10.219635550
    ref_ΩD = 0.055212384
    ref_ΩOD = 0.968016319
    ref_Ωtotal = 11.242864254
    @test wout.ΩI ≈ ref_ΩI
    @test wout.ΩD ≈ ref_ΩD
    @test wout.ΩOD ≈ ref_ΩOD
    @test wout.Ωtotal ≈ ref_Ωtotal
end

@testitem "read wout fortran stars" begin
    woutdir = joinpath(@__DIR__, "wout_testfiles")
    wout = read_wout(joinpath(woutdir, "stars.wout"))

    ref_centers = [
        [-0.866253, 1.973841, 1.973841],
        [-0.866253, 0.866253, 0.866253],
        [-99.973841, 1.973841, 0.866253],
        [NaN, 0.866253, 1.973841],
    ]
    @test wout.centers[1:3] ≈ ref_centers[1:3]
    @test isnan(wout.centers[end][1])
    @test wout.centers[end][2:end] ≈ ref_centers[end][2:end]
end
