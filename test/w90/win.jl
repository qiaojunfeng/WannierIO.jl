@testitem "read win" begin
    using LazyArtifacts
    win = read_win(artifact"Si2_valence/Si2_valence.win")

    # how to generate reference data:
    WRITE_TOML = false
    WRITE_TOML && write_win("/tmp/Si2_valence.win.toml"; toml=true, win...)

    test_data = read_win(artifact"Si2_valence/outputs/Si2_valence.win.toml")
    @test win == test_data
end

@testitem "read/write win" begin
    using LazyArtifacts
    win = read_win(artifact"Si2_valence/Si2_valence.win")

    tmpfile = tempname(; cleanup=true)
    write_win(tmpfile; win...)
    win2 = read_win(tmpfile)
    # compare without order
    @test Dict(pairs(win)) == Dict(pairs(win2))
end

@testitem "read/write win toml" begin
    using LazyArtifacts
    toml_path = artifact"Si2_valence/outputs/Si2_valence.win.toml"
    win = read_win(toml_path)

    tmpfile = tempname(; cleanup=true)
    write_win(tmpfile, WannierIO.Wannier90Toml(); win...)
    win2 = read_win(tmpfile)
    # compare without order
    @test Dict(pairs(win)) == Dict(pairs(win2))
end

@testitem "read win: special cases" begin
    windir = joinpath(@__DIR__, "win_testfiles")
    for win in readdir(windir)
        if endswith(win, ".win")
            # If it doesn't throw exceptions, it's good enough
            @test read_win(joinpath(windir, win)) isa NamedTuple
        end
    end
end

@testitem "read win: unknown blocks" begin
    windir = joinpath(@__DIR__, "win_testfiles")
    win = read_win(joinpath(windir, "unknown_blocks.win"))

    @test win.unknown_a == ["A1", "A2"]
    @test win.unknown_b == ["B1 B2"]
end

@testitem "read win: atoms_cart bohr" begin
    windir = joinpath(@__DIR__, "win_testfiles")
    ref_win = read_win(joinpath(windir, "unknown_blocks.win"))
    win = read_win(joinpath(windir, "atoms_cart_bohr.win"))

    for ((atom1, pos1), (atom2, pos2)) in zip(win.atoms_frac, ref_win.atoms_frac)
        @test atom1 == atom2
        @test pos1 ≈ pos2
    end
end

@testitem "read win: exclude_bands" begin
    windir = joinpath(@__DIR__, "win_testfiles")
    win = read_win(joinpath(windir, "exclude_bands.win"))
    @test win.exclude_bands == [1, 3, 4, 5, 6]
end

@testitem "read win: explicit_kpath" begin
    using LazyArtifacts
    toml_path = artifact"GaAs/GaAs.win"
    win = read_win(toml_path)

    @test win.explicit_kpath_labels == [
        :L => [0.5, 0.5, 0.5],
        :G => [0.0, 0.0, 0.0],
        :X => [0.5, 0.0, 0.5],
        :X2 => [0.5, -0.5, 0.0],
        :K => [0.375, -0.375, 0.0],
    ]
    @test length(win.explicit_kpath) == 214
    @test win.explicit_kpath[[1, end-1]] == [
        [0.5, 0.5, 0.5],
        [0.012097, -0.012097, 0.0],
    ]

    tmpfile = tempname(; cleanup=true)
    write_win(tmpfile; win...)
    win2 = read_win(tmpfile)
    # compare without order
    @test Dict(pairs(win)) == Dict(pairs(win2))
end

@testitem "write win: OrderedDict" begin
    using LazyArtifacts, OrderedCollections
    win = read_win(artifact"GaAs/GaAs.win")

    # convert to OrderedDict
    win = OrderedDict(pairs(win))
    tmpfile = tempname(; cleanup=true)
    write_win(tmpfile, win)
    win2 = read_win(tmpfile)
    # compare without order
    @test Dict(pairs(win)) == Dict(pairs(win2))
end
