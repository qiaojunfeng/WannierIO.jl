@testitem "_parse_wout_lattice" begin
    file = joinpath(@__DIR__, "wout_testfiles", "lattice.txt")
    lines = readlines(file)
    lattice = WannierIO._parse_wout_lattice(lines)
    ref_lattice = [
        0.1 1.1 2.1
        0.2 0.2 0.2
        0.3 0.3 0.3
    ]
    @test lattice ≈ ref_lattice
end

@testitem "_parse_wout_recip_lattice" begin
    file = joinpath(@__DIR__, "wout_testfiles", "recip_lattice.txt")
    lines = readlines(file)
    recip_lattice = WannierIO._parse_wout_recip_lattice(lines)
    ref_recip_lattice = [
        0.1 1.1 2.1
        0.2 0.2 0.2
        0.3 0.3 0.3
    ]
    @test recip_lattice ≈ ref_recip_lattice
end

@testitem "_parse_wout_atoms" begin
    file = joinpath(@__DIR__, "wout_testfiles", "atoms.txt")
    lines = readlines(file)
    atom_labels, atom_positions = WannierIO._parse_wout_atoms(lines)
    ref_atom_labels = ["Si", "Si"]
    ref_atom_positions = [[0.00000, 0.00000, 0.00000], [0.25000, 0.25000, 0.25000]]
    @test atom_labels == ref_atom_labels
    @test atom_positions ≈ ref_atom_positions
end

@testitem "_parse_wout_disentangle" begin
    file = joinpath(@__DIR__, "wout_testfiles", "disentangle.txt")
    lines = readlines(file)
    results = WannierIO._parse_wout_disentangle(lines)
    @test results.iter == [1, 341, 342]
    @test results.ΩI_previous ≈ [25.38943399, 16.22884440, 16.22884440]
    @test results.ΩI_current ≈ [21.32896063, 16.22884440, 16.22884440]
    @test results.ΔΩI ≈ [1.904e-01, -1.883e-10, -1.799e-10]
end

@testitem "_parse_wout_wf_center_spread" begin
    file = joinpath(@__DIR__, "wout_testfiles", "wf_center_spread.txt")
    lines = readlines(file)
    centers, spreads, sum_centers, sum_spreads = WannierIO._parse_wout_wf_center_spread(
        lines
    )
    ref_centers = [[-0.000005, 0.000021, 0.000023], [0.000013, -0.000054, 0.000016]]
    ref_spreads = [2.56218734, 3.19493515]
    ref_sum_centers = [5.430528, 5.430529, 5.430534]
    ref_sum_spreads = 24.29446759
    @test centers ≈ ref_centers
    @test spreads ≈ ref_spreads
    @test sum_centers ≈ ref_sum_centers
    @test sum_spreads ≈ ref_sum_spreads
end

@testitem "_parse_wout_wannierize" begin
    file = joinpath(@__DIR__, "wout_testfiles", "wannierize.txt")
    lines = readlines(file)
    results = WannierIO._parse_wout_wannierize(lines)
    @test results.iter == [0, 1, 45]
    @test results.centers ≈ [
        [[-0.000005, 0.000021, 0.000023]],
        [[-0.000005, 0.000020, 0.000022]],
        [[0.000001, 0.000006, 0.000006]],
    ]
    @test results.spreads ≈ [[2.56218734], [2.46316318], [1.95373328]]
    @test results.sum_centers ≈ [
        [5.430528, 5.430529, 5.430534],
        [5.430528, 5.430529, 5.430534],
        [5.430528, 5.430529, 5.430534],
    ]
    @test results.sum_spreads ≈ [24.29446759, 24.07813875, 23.58343588]
    @test results.ΩD ≈ [0.2135529, 0.2218113, 0.2612087]
    @test results.ΩOD ≈ [7.8520707, 7.6274835, 7.0933833]
    @test results.Ωtotal ≈ [24.2944680, 24.0781392, 23.5834363]
end

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
    @test wout.kgrid == [6, 6, 6]
    @test wout.sum_centers ≈ [-0.000001, -0.000000, 0.000000]
    @test wout.sum_spreads ≈ 7.71671562

    ref_ΩI = 7.153329184
    ref_ΩD = 0.000000000
    ref_ΩOD = 0.563386215
    ref_Ωtotal = 7.716715399
    @test wout.ΩI ≈ ref_ΩI
    @test wout.ΩD ≈ ref_ΩD
    @test wout.ΩOD ≈ ref_ΩOD
    @test wout.Ωtotal ≈ ref_Ωtotal

    ref_phases = [
        0.996157 + 0.087588im,
        0.996157 + 0.087588im,
        0.996157 + 0.087588im,
        0.998869 + 0.047543im,
    ]
    ref_imre = [4.566451, 4.566481, 4.566335, 2.154381]
    @test wout.phase_factors ≈ ref_phases
    @test wout.im_re_ratios ≈ ref_imre
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

@testitem "read wout iterations" begin
    using LazyArtifacts
    wout = read_wout(artifact"Fe_soc/reference/Fe.wout"; iterations=true)

    @test wout.iterations.disentangle.iter == 1:7000
    @test wout.iterations.disentangle.ΩI_previous[[1, 7000]] ≈ [16.06412287, 10.21963555]
    @test wout.iterations.disentangle.ΩI_current[[1, 7000]] ≈ [15.60506153, 10.21963555]
    @test wout.iterations.disentangle.ΔΩI[[1, 7000]] ≈ [2.942e-02, 3.224e-10]

    @test wout.iterations.wannierize.iter == 0:2083
    @test wout.iterations.wannierize.centers[[1, 2084]] ≈ [
        [
            [0.032917, -0.046451, -0.535906],
            [0.024248, 0.103848, 0.410673],
            [0.712774, -0.046012, -0.005995],
            [-0.833291, -0.010927, 0.004327],
            [0.008355, 0.067909, -0.010136],
            [0.019012, -0.111345, 0.041994],
            [-0.000064, 0.000132, -0.000702],
            [0.000076, -0.000232, 0.000050],
            [-0.003011, -0.000129, 0.000164],
            [0.002979, 0.000022, 0.000530],
            [0.000151, -0.000186, 0.000358],
            [-0.000021, 0.000357, 0.000260],
            [0.000171, -0.000092, 0.000177],
            [0.000002, 0.000041, -0.000426],
            [-0.002417, 0.000028, -0.000420],
            [0.002224, 0.000170, 0.000076],
        ],
        [
            [-0.001227, 0.761822, -1.210363],
            [-0.001132, 0.761153, 1.211403],
            [1.434030, 0.674118, -0.223106],
            [-1.436981, -0.000412, 0.958820],
            [-0.002088, -1.435412, 0.479468],
            [-0.001066, -0.761483, -1.210896],
            [0.000330, -0.000261, 0.000357],
            [0.000085, -0.019916, -0.010210],
            [-0.000258, -0.000766, 0.000475],
            [0.000067, -0.000720, -0.001011],
            [0.000082, 0.020905, 0.010305],
            [0.000045, 0.019902, -0.010005],
            [-0.000815, 0.000310, -0.000449],
            [0.000117, -0.020922, 0.010472],
            [-0.000259, 0.000986, 0.001254],
            [0.000055, 0.000493, -0.000470],
        ],
    ]
    @test wout.iterations.wannierize.spreads[[1, 2084]] ≈ [
        [
            6.91546411,
            7.13787642,
            9.87200446,
            9.48728884,
            5.67519982,
            5.93428114,
            0.60839403,
            0.58059954,
            0.64462091,
            0.61126568,
            0.63125454,
            0.62139130,
            0.61309583,
            0.59680091,
            0.64713955,
            0.61960630,
        ],
        [
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
            0.44007660,
            0.45389206,
            0.48648695,
            0.44062982,
        ],
    ]
    @test wout.iterations.wannierize.sum_centers[[1, 2084]] ≈
        [[-0.035896, -0.042867, -0.094975], [-0.009014, -0.000203, 0.006043]]
    @test wout.iterations.wannierize.sum_spreads[[1, 2084]] ≈ [51.19628337, 11.24286425]
    @test wout.iterations.wannierize.ΩD[[1, 2084]] ≈ [22.8931526, 0.0552124]
    @test wout.iterations.wannierize.ΩOD[[1, 2084]] ≈ [18.0834952, 0.9680163]
    @test wout.iterations.wannierize.Ωtotal[[1, 2084]] ≈ [51.1962834, 11.2428643]

    # test that parsing the "Final State" is correct
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
