using WannierIO: Vec3

@testset "read qe xml" begin
    qe = WannierIO.read_qe_xml(joinpath(FIXTURE_PATH, "qe/si2.xml"))

    lattice = [0.0 2.715265 2.715265; 2.715265 0.0 2.715265; 2.715265 2.715265 0.0]
    @test qe.lattice ≈ lattice

    atom_positions = [Vec3(0.0, 0.0, 0.0), Vec3(0.25, 0.25, 0.25)]
    @test qe.atom_positions ≈ atom_positions

    atom_labels = ["Si", "Si"]
    @test qe.atom_labels == atom_labels

    recip_lattice = [
        -1.1570114348285685 1.1570114348285685 1.1570114348285685
        1.1570114348285685 -1.1570114348285685 1.1570114348285685
        1.1570114348285685 1.1570114348285685 -1.1570114348285685
    ]
    @test qe.recip_lattice ≈ recip_lattice

    kpoints = Vec3[
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.5],
        [0.0, 0.5, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.0],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
        [0.5, 0.5, 0.5],
    ]
    @test qe.kpoints ≈ kpoints

    eigenvalues = [
        [-5.798515544174599, 6.202325333825356, 6.202325333825502, 6.20232533387998],
        [-3.445835420133952, -0.8272627456868693, 4.993059918379725, 4.993059918395274],
        [-3.445835420133764, -0.8272627456869411, 4.993059918286839, 4.993059918321653],
        [-1.6353082316074206, -1.6353082316054404, 3.3159306458387663, 3.315930645895712],
        [-3.4458354201337613, -0.8272627456869834, 4.993059918283482, 4.993059918316889],
        [-1.6353082316050855, -1.6353082316009113, 3.3159306458027276, 3.3159306458796163],
        [-1.6353082316059693, -1.6353082315998204, 3.315930645843575, 3.3159306458994777],
        [-3.445835420133952, -0.8272627456862578, 4.993059918290434, 4.993059918308719],
    ]
    @test qe.eigenvalues ≈ eigenvalues

    fermi_energy = 7.1028327867214625
    @test qe.fermi_energy ≈ fermi_energy
end

@testset "read qe xml spin-polarized" begin
    qe = WannierIO.read_qe_xml(joinpath(FIXTURE_PATH, "qe/cri3.xml"))

    lattice = [
        6.8171434485254725 -3.4085717242627362 0.0
        0.0 5.903819407666132 0.0
        0.0 0.0 20.078373841305446
    ]
    @test qe.lattice ≈ lattice

    atom_positions = Vec3[
        [0.33333333330000003, 0.6666666667000001, 0.0],
        [0.6666666667, 0.3333333332999999, 0.0],
        [0.0, 0.35410773599999995, 0.0769313053],
        [0.0, 0.645892264, 0.9230687247000002],
        [0.354107736, 0.35410773599999995, 0.9230687247000002],
        [0.645892264, 0.645892264, 0.0769313053],
        [0.354107736, 0.0, 0.0769313053],
        [0.645892264, 0.0, 0.9230687247000002],
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

    kpoints = Vec3[[0.0, 0.0, 0.0], [0.0, 0.5, 0.0]]
    @test qe.kpoints ≈ kpoints

    eigenvalues_up = [
        [
            -74.54574717483933,
            -74.54551689417254,
            -46.011120064567244,
            -46.010179445079764,
            -46.00264241832006,
            -46.00264241780246,
            -45.99945365983328,
            -45.999453659349406,
            -15.258137817774953,
            -14.434129479400292,
            -14.128321864637826,
            -14.12832185282107,
            -14.002018642981103,
            -14.00201863176426,
            -6.955103201948282,
            -6.377830749851818,
            -5.832064885668679,
            -5.832064877937429,
            -5.785438613929024,
            -5.785438609151065,
            -5.470966913010504,
            -4.985335808178225,
            -4.98533580745564,
            -4.916951525453973,
            -4.499572349388002,
            -4.499572342557247,
            -3.9947232080300425,
            -3.4100364649081154,
            -3.410036462682431,
            -3.3432369930835293,
            -3.212853657810935,
            -2.877148108167079,
            -2.8771481072714242,
            -2.8205377091155506,
            -2.7112546034220335,
            -2.7112546001204527,
            -2.227864302964817,
            -2.227864302681989,
            -0.9649879408667482,
            -0.9649879386921768,
            -0.5011821567286011,
            -0.5011821555965684,
        ],
        [
            -74.54567100208304,
            -74.54559247336806,
            -46.01080823523272,
            -46.01049522953251,
            -46.00284774500801,
            -46.001786606172075,
            -46.00030942778289,
            -45.99924502288391,
            -14.866260631660657,
            -14.60830302443693,
            -14.226086190297028,
            -14.21719317285224,
            -14.139950296193405,
            -14.047361593591829,
            -6.875063642828858,
            -6.188090873134216,
            -6.008584566449626,
            -5.7503990828336775,
            -5.497569331219267,
            -5.46514213633486,
            -5.352052023769881,
            -5.009062562591209,
            -4.813489491306931,
            -4.738629649860174,
            -4.503559190643073,
            -4.382613668875345,
            -4.060468433018202,
            -3.7011838264277976,
            -3.5240503142480706,
            -3.500776444834226,
            -3.2246759737161876,
            -3.1048638329192637,
            -2.987196470750001,
            -2.7296169076374253,
            -2.7030296916618033,
            -2.527174855492087,
            -2.50162144572179,
            -2.249445541011562,
            -0.8801725163288753,
            -0.8796264762539415,
            -0.8649898770312278,
            -0.7184951144043287,
        ],
    ]
    @test qe.eigenvalues_up ≈ eigenvalues_up

    eigenvalues_dn = [
        [
            -71.39697269829512,
            -71.39669042491292,
            -42.93568753391549,
            -42.93434357617464,
            -42.89798663678183,
            -42.897986636259695,
            -42.89381069743494,
            -42.89381069695789,
            -15.246164375300854,
            -14.470204933679115,
            -14.161594209549216,
            -14.161594197790915,
            -14.035765615475894,
            -14.035765604276474,
            -6.647596352761345,
            -6.1222167267291265,
            -5.658786285385981,
            -5.6587862774004485,
            -5.507676957206841,
            -5.507676952192285,
            -5.398831415094653,
            -4.686652758666972,
            -4.372689958861119,
            -4.372689956596873,
            -4.017309103202073,
            -4.017309093993718,
            -3.9512527770398354,
            -3.3644935561853435,
            -2.9839135141316007,
            -2.9839135101493204,
            -2.405543369382001,
            -2.4055433685404344,
            -0.2901608432110109,
            -0.29016084245832685,
            -0.2830193408701903,
            -0.21416631711254216,
            0.07715702774977284,
            0.07715702798788102,
            0.5365081337392321,
            0.5365081346986672,
            0.856487375430603,
            0.8564873762311315,
        ],
        [
            -71.3968747216541,
            -71.39677918611643,
            -42.93524317519616,
            -42.93479498786359,
            -42.89826946646218,
            -42.89688132082972,
            -42.894916264901674,
            -42.89352337405696,
            -14.864089135674183,
            -14.622783931130696,
            -14.259934393459213,
            -14.253368483029417,
            -14.173640487596431,
            -14.078307770947594,
            -6.617138407087747,
            -5.929492886946907,
            -5.929251953835399,
            -5.393642396197951,
            -5.262558098151388,
            -5.256162971712015,
            -4.914857535424112,
            -4.5896635668867605,
            -4.38320083135312,
            -4.199838339973026,
            -4.134745675917602,
            -3.9932448624185146,
            -3.7202712178030977,
            -3.4427700318108596,
            -3.394790639926803,
            -3.2359060490074647,
            -2.9989871386326685,
            -2.6206639901118365,
            -0.4837327947227868,
            -0.44883507577724074,
            -0.14922557722736784,
            0.013746024341250367,
            0.04363885877401906,
            0.15203798491316595,
            0.477526367533286,
            0.5833090750983857,
            0.7162584520246695,
            0.7310051235665301,
        ],
    ]
    @test qe.eigenvalues_dn ≈ eigenvalues_dn

    fermi_energy = -1.3722789296611695
    @test qe.fermi_energy ≈ fermi_energy
end
