@testitem "read qe xml" begin
    using WannierIO: Vec3
    using LazyArtifacts
    qe = WannierIO.read_qe_xml(artifact"Si2/reference/qe_bands.xml")

    lattice = [0.0 2.715265 2.715265; 2.715265 0.0 2.715265; 2.715265 2.715265 0.0]
    @test qe.lattice ≈ lattice

    atom_positions = Vec3[[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]]
    @test qe.atom_positions ≈ atom_positions

    atom_labels = ["Si", "Si"]
    @test qe.atom_labels == atom_labels

    recip_lattice = [
        -1.1570114348285685 1.1570114348285685 1.1570114348285685
        1.1570114348285685 -1.1570114348285685 1.1570114348285685
        1.1570114348285685 1.1570114348285685 -1.1570114348285685
    ]
    @test qe.recip_lattice ≈ recip_lattice

    @test length(qe.kpoints) == 511
    @test qe.kpoints[1] ≈ [0.0, 0.0, 0.0]
    @test qe.kpoints[end] ≈ [0.5, 0.0, 0.5]

    @test length(qe.eigenvalues) == 511
    eigenvalues1 = [
        -5.826225550687528,
        6.165602316381265,
        6.165602316381265,
        6.165602316381303,
        8.683797578503912,
        8.683797578503947,
        8.683797578503983,
        9.477127583014145,
        13.86461804570024,
        13.864618045701016,
        13.87950848504097,
        17.27314636025589,
        17.273146360255925,
        17.273146360257257,
        21.254697588478756,
        29.07002686142792,
    ]
    eigenvalues511 = [
        -1.6666776531834973,
        -1.6666776531834973,
        3.2884932258593533,
        3.2884932258593533,
        6.758325832996588,
        6.758325832996588,
        16.220036706342146,
        16.220036706342146,
        17.107712425166522,
        17.10771242516657,
        18.670405795535505,
        18.67040579553565,
        19.10855942754898,
        19.10855942754903,
        24.841843066432848,
        24.841843066563083,
    ]
    @test qe.eigenvalues[1] ≈ eigenvalues1
    @test qe.eigenvalues[end] ≈ eigenvalues511

    @test qe.n_electrons ≈ 8.0
    @test qe.fermi_energy ≈ 6.528341904366175
    @test qe.alat ≈ 3.8399645884368714
end

@testitem "read qe xml spin-polarized" begin
    using WannierIO: Vec3
    using LazyArtifacts
    qe = WannierIO.read_qe_xml(artifact"CrI3/reference/qe_bands.xml")

    lattice = [
        6.8171434485254725 -3.4085717242627362 0.0
        0.0 5.903819407666132 0.0
        0.0 0.0 20.078373841305446
    ]
    @test qe.lattice ≈ lattice

    atom_positions = Vec3[
        [0.33333333330000003, 0.6666666667000001, 0.0],
        [0.6666666667, 0.3333333332999999, 0.0],
        [-2.7755575615628914e-17, 0.35410773599999995, 0.0769313053],
        [0.0, 0.645892264, -0.07693127530000003],
        [0.354107736, 0.35410773599999995, -0.07693127530000003],
        [0.645892264, 0.645892264, 0.0769313053],
        [0.354107736, 0.0, 0.0769313053],
        [0.645892264, 0.0, -0.07693127530000003],
    ]
    @test qe.atom_positions ≈ atom_positions

    atom_labels = ["Cr", "Cr", "I", "I", "I", "I", "I", "I"]
    @test qe.atom_labels == atom_labels

    recip_lattice = [
        0.921674210704576 0.0 0.0
        0.532128853655385 1.0642577073107704 0.0
        0.0 0.0 0.31293297738354436
    ]
    @test qe.recip_lattice ≈ recip_lattice

    @test length(qe.kpoints) == 274
    kpoint2 = Vec3(0.0, 0.0049999999999999975, 0.0)
    @test qe.kpoints[2] ≈ kpoint2

    @test length(qe.eigenvalues_up) == 274
    eigenvalues_up2 = [-77.99192823029188, -77.99169805183234, -49.45736655071318]
    @test qe.eigenvalues_up[2][1:3] ≈ eigenvalues_up2

    @test length(qe.eigenvalues_dn) == 274
    eigenvalues_dn2 = [-74.843061971795, -74.84277910814951, -46.38172392618895]
    @test qe.eigenvalues_dn[2][1:3] ≈ eigenvalues_dn2

    @test qe.fermi_energy ≈ -4.819375066024118
end
