using YAML

@testset "read win" begin
    test_data = YAML.load_file(String(@__DIR__) * "/test_data/win.yaml")

    win = read_win(joinpath(FIXTURE_PATH, "si2.win"))

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
    #     "kpoint_path" => win.kpoint_path,
    # )
    # YAML.write_file(String(@__DIR__) * "/test_data/win.yaml", yaml_dict)

    function test_kpoint_path()
        test_kpath = test_data["kpoint_path"]
        for (i, path) in enumerate(test_kpath)
            win_path = win.kpoint_path[i]
            k1, k2 = path
            # Dict to pair
            k1 = first(k1)
            k2 = first(k2)

            (Symbol(k1.first) != win_path[1].first) && return false
            (k1.second ≉ win_path[1].second) && return false
            (Symbol(k2.first) != win_path[2].first) && return false
            (k2.second ≉ win_path[2].second) && return false
        end
        return true
    end
    @test test_kpoint_path()

    @test test_data["num_wann"] == win.num_wann
    @test test_data["num_bands"] == win.num_bands
    @test test_data["unit_cell"] ≈ unit_cell
    @test test_data["atoms_frac"] ≈ atoms_frac
    @test test_data["atom_labels"] == win.atom_labels
    @test test_data["mp_grid"] == win.mp_grid
    @test test_data["kpoints"] ≈ kpoints
    @test ismissing(win.dis_froz_min)
    @test ismissing(win.dis_win_min)
    @test test_data["dis_froz_max"] == "missing"
    @test test_data["dis_win_max"] == "missing"
end
