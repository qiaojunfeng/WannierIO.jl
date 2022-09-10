
@testset "read/write unk" begin
    ik, Ψ = read_unk(joinpath(FIXTURE_PATH, "silicon/UNK00001.1"))

    tmpfile = tempname(; cleanup=true)
    write_unk(tmpfile, ik, Ψ)
    ik2, Ψ2 = read_unk(tmpfile)

    @test ik ≈ ik2
    @test Ψ ≈ Ψ2
end
