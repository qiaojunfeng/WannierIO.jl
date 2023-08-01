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
    using WannierIO: Vec3
    using LazyArtifacts
    wout = read_wout(artifact"Fe/reference/Fe.wout")

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
    ref_atom_positions = Vec3[[0.00000, 0.00000, 0.00000]]
    ref_centers = Vec3[
        [0.258343, -0.648816, -1.435000],
        [-0.256199, 0.651431, 1.434999],
        [-1.177202, -0.785555, -0.000002],
        [1.177057, 0.785544, 0.000001],
        [-0.556737, 1.434443, -0.000006],
        [0.560726, -1.435175, 0.000007],
        [-0.000012, -0.000007, -0.000074],
        [-0.000039, 0.000388, 0.000024],
        [0.000193, 0.000181, 0.000073],
        [-0.000202, -0.000184, -0.000012],
        [0.000089, -0.000083, 0.000044],
        [-0.000068, 0.000065, -0.000008],
        [0.020630, -0.001504, -0.000001],
        [-0.020451, -0.001798, -0.000029],
        [-0.000180, 0.002411, -0.000005],
        [0.000124, 0.000352, -0.000011],
    ]
    ref_spreads = [
        1.23964450,
        1.23815613,
        1.23964003,
        1.23812433,
        1.26093796,
        1.26187450,
        0.43635868,
        0.39242894,
        0.45563745,
        0.41876201,
        0.45572233,
        0.41876485,
        0.47408733,
        0.47011636,
        0.40823901,
        0.39315249,
    ]

    @test wout.lattice ≈ ref_lattice
    @test wout.recip_lattice ≈ ref_recip_lattice
    @test wout.atom_labels == ref_atom_labels
    @test wout.atom_positions ≈ ref_atom_positions
    @test wout.centers ≈ ref_centers
    @test wout.spreads ≈ ref_spreads

    ref_ΩI = 10.545296062
    ref_ΩD = 0.064753479
    ref_ΩOD = 1.191596786
    ref_Ωtotal = 11.801646328
    @test wout.ΩI ≈ ref_ΩI
    @test wout.ΩD ≈ ref_ΩD
    @test wout.ΩOD ≈ ref_ΩOD
    @test wout.Ωtotal ≈ ref_Ωtotal
end
