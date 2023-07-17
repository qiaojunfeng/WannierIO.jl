@testitem "read win" begin
    using LazyArtifacts
    win = read_win(artifact"Si2_valence/Si2_valence.win")

    # how to generate reference data:
    WRITE_TOML = false
    WRITE_TOML && write_win("/tmp/Si2_valence.win.toml"; toml=true, win...)

    test_data = read_win(artifact"Si2_valence/reference/Si2_valence.win.toml")
    @test win == test_data
end

@testitem "read/write win" begin
    using LazyArtifacts
    win = read_win(artifact"Si2_valence/Si2_valence.win")

    tmpfile = tempname(; cleanup=true)
    write_win(tmpfile; win...)
    win2 = read_win(tmpfile)
    @test win == win2
end

@testitem "read/write win toml" begin
    using LazyArtifacts
    toml_path = artifact"Si2_valence/reference/Si2_valence.win.toml"
    win = read_win(toml_path)

    tmpfile = tempname(; cleanup=true)
    write_win(tmpfile; win...)
    win2 = read_win(tmpfile)
    @test win == win2
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
