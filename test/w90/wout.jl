using WannierIO: Vec3
@testset "read wout" begin
    wout = read_wout(joinpath(FIXTURE_PATH, "si2.wout"))

    ref_lattice =
        [
            0.000000 2.715265 2.715265
            2.715265 0.000000 2.715265
            2.715265 2.715265 0.000000
        ]'
    ref_atom_labels = ["Si", "Si"]
    ref_atom_positions = [Vec3(0.00000, 0.00000, 0.00000), Vec3(0.25000, 0.25000, 0.25000)]
    ref_centers = [
        Vec3(0.678817, -0.678817, -0.678816),
        Vec3(-0.678817, -0.678816, 0.678816),
        Vec3(-0.678816, 0.678816, -0.678816),
        Vec3(0.678816, 0.678816, 0.678817),
    ]
    ref_spreads = [
        1.02317854
        1.02317859
        1.02317884
        1.02317863
    ]

    @test wout.lattice ≈ ref_lattice
    @test wout.atom_labels == ref_atom_labels
    @test wout.atom_positions ≈ ref_atom_positions
    @test wout.centers ≈ ref_centers
    @test wout.spreads ≈ ref_spreads
end
