@testset "read wout" begin
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
end

# TODO: add test for non symmetric lattice, with disentanglement?
