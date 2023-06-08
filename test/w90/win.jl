@testset "read win" begin
    win = read_win(artifact"Si2_valence/si2.win")
    toml_path = artifact"Si2_valence/reference/si2.win.toml"

    # how to generate reference data:
    WRITE_TOML = false
    WRITE_TOML && write_win(toml_path, win; toml=true)

    test_data = read_win(toml_path)
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
