@testset "read chk" begin
    chk = read_chk(artifact"Si2_valence/reference/Si2_valence.chk.fmt")

    @test chk.n_wann == 4
    @test chk.n_bands == 4
    @test chk.n_kpts == 6^3
    @test chk.n_bvecs == 8
    @test length(chk.M) == 6^3
    @test length(chk.M[1]) == 8
    @test size(chk.M[1][1]) == (4, 4)
    ref_M12 = ComplexF64[
        0.9170561763002493-0.12073283760315105im 0.006084376539889745-0.0008009985221228121im 0.06241876168424332+0.02585491031881253im 0.006084587068619811-0.0008011618858289517im
        0.006084880213061389-0.0008013918486074993im 0.9170555712700101-0.1207327246826146im 0.06241866709702579+0.025854802566432057im 0.006084809499489799-0.000800887926406555im
        0.04089804920678168-0.005384004515382356im 0.040899090696488685-0.005384672958981166im 0.789621736965331+0.32707172642290266im 0.04089840122633419-0.005384316593360309im
        0.006084658804818147-0.0008009811562424997im 0.0060846064084139105-0.0008007905091687038im 0.0624189360729618+0.02585459811803456im 0.9170558895434384-0.1207329088207673im
    ]
    @test chk.M[1][2] ≈ ref_M12

    @test length(chk.Uml) == 216
    @test size(chk.Uml[1]) == (4, 4)
    ref_U2 = ComplexF64[
        -0.16235109810009157+0.47353334937462366im -0.3755805652160198+0.32735938468233533im -0.16235110057874652+0.4735333026272855im -0.16235112363576365+0.47353331917208363im
        0.16956618931090092-0.23235535713863867im -0.7928314959724974+0.35098901875473565im 0.16956603006099888-0.23235524733459334im 0.16956607268060825-0.232355268156885im
        -0.3123360136432703+0.025184486826638694im -3.782746554514418e-8+1.6862208020401383e-8im -0.10111009556353681-0.6127313622218364im 0.41344612298635436+0.587546837211439im
        -0.6776552858574172-0.33054780770943226im -8.7307707370288e-8-6.01196661668837e-8im 0.30384128248067205+0.43437798438065095im 0.3738141419013159-0.10383009885606777im
    ]
    @test chk.Uml[2] ≈ ref_U2

    @test chk.Udis == Matrix{ComplexF64}[]
    @test chk.checkpoint == "postwann"
    @test chk.dis_bands == BitVector[]
    @test chk.n_dis == Int64[]
    @test chk.exclude_bands == Int64[]
    @test chk.n_exclude_bands == 0
    @test chk.have_disentangled == false
    @test chk.header == "written on 15Jun2023 at 10:39:45"

    win = read_win(artifact"Si2_valence/Si2_valence.win")
    @test chk.kgrid == win.mp_grid
    @test chk.kpoints == win.kpoints

    # make sure we read the lattice as column-major
    @test chk.lattice ≈ win.unit_cell_cart

    @test chk.recip_lattice ≈ WannierIO.get_recip_lattice(chk.lattice)

    ref_r = WannierIO.Vec3[
        [0.6788160683908873, -0.678816205763796, -0.6788162184647731],
        [-0.6788164049684758, -0.6788163321814957, 0.6788162805858188],
        [-0.6788162635506181, 0.6788163021769287, -0.678816181780838],
        [0.6788160891131698, 0.6788160528597309, 0.6788162387951322],
    ]
    @test chk.r ≈ ref_r
    @test chk.ΩI == -1.0
    ref_ω = [1.9291789239559392, 1.9291789959111603, 1.9291788335828866, 1.9291788667547383]
    @test chk.ω ≈ ref_ω
end

@testset "read/write chk" begin
    chk = read_chk(artifact"Si2_valence/reference/Si2_valence.chk.fmt")

    tmpfile = tempname(; cleanup=true)
    write_chk(tmpfile, chk)
    chk2 = read_chk(tmpfile)
    @test chk ≈ chk2
end

@testset "read/write chk binary" begin
    chk = read_chk(artifact"Si2_valence/reference/Si2_valence.chk.fmt")
    chk1 = read_chk(artifact"Si2_valence/reference/binary/Si2_valence.chk")
    @test chk ≈ chk1

    tmpfile = tempname(; cleanup=true)
    write_chk(tmpfile, chk; binary=true)
    chk2 = read_chk(tmpfile)
    @test chk ≈ chk2
end

@testset "read chk disentanglement" begin
    chk = read_chk(artifact"Si2/reference/Si2.chk.fmt")

    @test chk.n_wann == 8
    @test chk.n_bands == 16
    @test chk.n_kpts == 9^3
    @test chk.n_bvecs == 8
    @test length(chk.M) == 9^3
    @test length(chk.M[1]) == 8
    @test size(chk.M[1][1]) == (8, 8)
    @test chk.dis_bands == [trues(16) for _ in 1:(9^3)]
    @test chk.n_dis == [16 for _ in 1:(9^3)]
    @test chk.have_disentangled == true
end

@testset "read/write chk disentanglement" begin
    chk = read_chk(artifact"Si2/reference/Si2.chk.fmt")

    tmpfile = tempname(; cleanup=true)
    write_chk(tmpfile, chk)
    chk2 = read_chk(tmpfile)
    @test chk ≈ chk2
end

@testset "read/write chk disentanglement binary" begin
    chk = read_chk(artifact"Si2/reference/Si2.chk.fmt")
    chk1 = read_chk(artifact"Si2/reference/Si2.chk")
    @test chk ≈ chk1

    tmpfile = tempname(; cleanup=true)
    write_chk(tmpfile, chk; binary=true)
    chk2 = read_chk(tmpfile)
    @test chk ≈ chk2
end

@testset "get_Udis" begin
    chk = read_chk(artifact"Si2/reference/Si2.chk")
    Udis = get_Udis(chk)
    Udis_ref = read_amn(artifact"Si2/reference/Si2.chk_Udis.amn")
    @test Udis ≈ Udis_ref
end

@testset "get_U" begin
    chk = read_chk(artifact"Si2/reference/Si2.chk")
    U = get_U(chk)
    U_ref = read_amn(artifact"Si2/reference/Si2.chk_U.amn")
    @test U ≈ U_ref
end

# TODO: test chk recip lattice column-major
