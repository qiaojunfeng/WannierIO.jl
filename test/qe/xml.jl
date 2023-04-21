@testset "read qe xml" begin
    qe = WannierIO.read_qe_xml(joinpath(FIXTURE_PATH, "qe/si2.xml"))

    lattice = [0.0 2.715265 2.715265; 2.715265 0.0 2.715265; 2.715265 2.715265 0.0]
    @test qe.lattice ≈ lattice

    atom_positions = [0.0 0.25; 0.0 0.25; 0.0 0.25]
    @test qe.atom_positions ≈ atom_positions

    atom_labels = ["Si", "Si"]
    @test qe.atom_labels == atom_labels

    recip_lattice = [
        -1.1570114348285685 1.1570114348285685 1.1570114348285685
        1.1570114348285685 -1.1570114348285685 1.1570114348285685
        1.1570114348285685 1.1570114348285685 -1.1570114348285685
    ]
    @test qe.recip_lattice ≈ recip_lattice

    kpoints = [
        0.0 0.0 0.0 0.0 0.5 0.5 0.5 0.5
        0.0 0.0 0.5 0.5 0.0 0.0 0.5 0.5
        0.0 0.5 0.0 0.5 0.0 0.5 0.0 0.5
    ]
    @test qe.kpoints ≈ kpoints

    E = [
        -5.798515544174599 -3.445835420133952 -3.445835420133764 -1.6353082316074206 -3.4458354201337613 -1.6353082316050855 -1.6353082316059693 -3.445835420133952
        6.202325333825356 -0.8272627456868693 -0.8272627456869411 -1.6353082316054404 -0.8272627456869834 -1.6353082316009113 -1.6353082315998204 -0.8272627456862578
        6.202325333825502 4.993059918379725 4.993059918286839 3.3159306458387663 4.993059918283482 3.3159306458027276 3.315930645843575 4.993059918290434
        6.20232533387998 4.993059918395274 4.993059918321653 3.315930645895712 4.993059918316889 3.3159306458796163 3.3159306458994777 4.993059918308719
    ]
    @test qe.E ≈ E

    fermi_energy = 7.1028327867214625
    @test qe.fermi_energy ≈ fermi_energy
end