
@testset "read wout" begin
    wout = read_wout(joinpath(FIXTURE_PATH, "valence/band/silicon.wout"))

    ref_lattice =
        [
            -2.698804 0.000000 2.698804
            0.000000 2.698804 2.698804
            -2.698804 2.698804 0.000000
        ]'
    ref_atom_labels = ["Si", "Si"]
    ref_atom_positions = [
        -0.25000 0.75000 -0.25000
        0.00000 0.00000 0.00000
    ]'
    ref_centers =
        [
            -0.659352 0.658238 -0.680969
            0.669283 0.695828 0.666806
            0.682490 -0.683846 -0.683726
            -0.701673 -0.656575 0.703751
        ]'
    ref_spreads = [
        2.39492617
        2.19372718
        1.83863803
        1.88512458
    ]

    @test wout.lattice ≈ ref_lattice
    @test wout.atom_labels == ref_atom_labels
    @test wout.atom_positions ≈ ref_atom_positions
    @test wout.centers ≈ ref_centers
    @test wout.spreads ≈ ref_spreads
end
