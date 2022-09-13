function Base.isapprox(a::WannierIO.Chk, b::WannierIO.Chk)
    for f in fieldnames(typeof(a))
        if getfield(a, f) isa String
            getfield(a, f) == getfield(b, f) || return false
        else
            getfield(a, f) ≈ getfield(b, f) || return false
        end
    end
    return true
end

@testset "read chk" begin
    chk = read_chk(joinpath(FIXTURE_PATH, "formatted/si2.chk.fmt"))

    @test chk.n_wann == 4
    @test chk.n_bands == 4

    # make sure we read the lattice as column-major
    win = read_win(joinpath(FIXTURE_PATH, "si2.win"))
    @test chk.lattice ≈ win.unit_cell
end

@testset "read/write chk" begin
    chk = read_chk(joinpath(FIXTURE_PATH, "formatted/si2.chk.fmt"))

    tmpfile = tempname(; cleanup=true)
    write_chk(tmpfile, chk)
    chk2 = read_chk(tmpfile)

    @test chk ≈ chk2
end

@testset "read/write chk binary" begin
    chk = read_chk(joinpath(FIXTURE_PATH, "formatted/si2.chk.fmt"))
    chk1 = read_chk(joinpath(FIXTURE_PATH, "unformatted/si2.chk"))

    @test chk ≈ chk1

    tmpfile = tempname(; cleanup=false)#cleanup=true)
    write_chk(tmpfile, chk; binary=true)
    chk2 = read_chk(tmpfile)

    @test chk ≈ chk2
end
