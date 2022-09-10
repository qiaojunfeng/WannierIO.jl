using YAML

mat2vec(A::AbstractMatrix) = [v for v in eachcol(A)]

@testset "read win" begin
    test_data = YAML.load_file(String(@__DIR__) * "/test_data/win.yaml")

    win = read_win(joinpath(FIXTURE_PATH, "silicon/silicon.win"))

    # Convert type so YAML can write it.
    unit_cell = mat2vec(win.unit_cell)
    atoms_frac = mat2vec(win.atoms_frac)
    kpoints = mat2vec(win.kpoints)

    # yaml_dict = Dict(
    #     "num_bands" => win.num_bands,
    #     "num_wann" => win.num_wann,
    #     "atoms_frac" => atoms_frac,
    #     "atom_labels" => win.atom_labels,
    #     "unit_cell" => unit_cell,
    #     "mp_grid" => win.mp_grid,
    #     "kpoints" => kpoints,
    #     "dis_froz_max" => win.dis_froz_max,
    #     "dis_win_max" => win.dis_win_max,
    #     "kpoint_path" => Dict(
    #         "basis" => win.kpoint_path.basis,
    #         "paths" => win.kpoint_path.paths,
    #         "points" => win.kpoint_path.points,
    #         "setting" => win.kpoint_path.setting,
    #     ),
    # )
    # YAML.write_file(String(@__DIR__) * "/test_data/win.yaml", yaml_dict)

    test_kpath = test_data["kpoint_path"]
    win_kpath = win.kpoint_path
    t_points = Dict(Symbol(k) => v for (k, v) in test_kpath["points"])
    @test t_points == win_kpath.points
    t_paths = [[Symbol(i) for i in l] for l in test_kpath["paths"]]
    @test t_paths == win_kpath.paths
    @test test_kpath["basis"] == win_kpath.basis
    @test Symbol(test_kpath["setting"]) == Symbol(win_kpath.setting)

    @test test_data["num_wann"] == win.num_wann
    @test test_data["num_bands"] == win.num_bands
    @test test_data["unit_cell"] ≈ unit_cell
    @test test_data["atoms_frac"] ≈ atoms_frac
    @test test_data["atom_labels"] == win.atom_labels
    @test test_data["mp_grid"] == win.mp_grid
    @test test_data["kpoints"] ≈ kpoints
    @test ismissing(win.dis_froz_min)
    @test ismissing(win.dis_win_min)
    @test test_data["dis_froz_max"] == win.dis_froz_max
    @test test_data["dis_win_max"] == win.dis_win_max
end
