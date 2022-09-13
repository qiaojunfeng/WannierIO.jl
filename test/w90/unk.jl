
@testset "read/write unk" begin
    ik, Ψ = read_unk(joinpath(FIXTURE_PATH, "formatted/UNK00001.1"))

    tmpfile = tempname(; cleanup=true)
    write_unk(tmpfile, ik, Ψ)
    ik2, Ψ2 = read_unk(tmpfile)

    @test ik ≈ ik2
    @test Ψ ≈ Ψ2
end

@testset "read/write unk binary" begin
    ik, Ψ = read_unk(joinpath(FIXTURE_PATH, "formatted/UNK00001.1"))
    ik1, Ψ1 = read_unk(joinpath(FIXTURE_PATH, "unformatted/UNK00001.1"))
    @test ik == ik1
    @test Ψ ≈ Ψ1

    tmpfile = tempname(; cleanup=true)
    write_unk(tmpfile, ik, Ψ; binary=true)
    ik2, Ψ2 = read_unk(tmpfile)

    @test ik == ik2
    @test Ψ ≈ Ψ2
end
