@testset "read win" begin
    win = read_win(artifact"Si2_valence/si2.win")

    # how to generate reference data:
    WRITE_TOML = false
    WRITE_TOML && WannierIO._write_win_toml("/tmp/si2.win.toml"; win...)

    test_data = WannierIO._read_win_toml(artifact"Si2_valence/reference/si2.win.toml")

    @test win == test_data
end

@testset "read/write win" begin
    win = read_win(artifact"Si2_valence/si2.win")

    tmpfile = tempname(; cleanup=true)
    write_win(tmpfile; win...)
    win2 = read_win(tmpfile)
    @test win == win2
end

@testset "read/write win toml" begin
    toml_path = artifact"Si2_valence/reference/si2.win.toml"
    win = WannierIO._read_win_toml(toml_path)

    tmpfile = tempname(; cleanup=true)
    write_win(tmpfile; win...)
    win2 = read_win(tmpfile)
    @test win == win2
end

@testset "read win: special cases" begin
    windir = joinpath(@__DIR__, "win/")
    for win in readdir(windir)
        if endswith(win, ".win")
            # If it doesn't throw exceptions, it's good enough
            @test read_win(joinpath(windir, win)) isa NamedTuple
        end
    end
end
