
@testset "read/write amn" begin
    A = read_amn(joinpath(FIXTURE_PATH, "formatted/si2.amn"))

    tmpfile = tempname(; cleanup=true)

    write_amn(tmpfile, A)

    A2 = read_amn(tmpfile)

    @test A ≈ A2
end
