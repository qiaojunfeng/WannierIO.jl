@testset "read/write amn" begin
    A = read_amn(joinpath(FIXTURE_PATH, "formatted/si2.amn"))

    tmpfile = tempname(; cleanup=true)
    write_amn(tmpfile, A)
    A2 = read_amn(tmpfile)
    @test A ≈ A2
end

@testset "read/write amn binary" begin
    A = read_amn(joinpath(FIXTURE_PATH, "formatted/si2.amn"))
    A1 = read_amn(joinpath(FIXTURE_PATH, "unformatted/si2.amn"))
    @test A ≈ A1

    tmpfile = tempname(; cleanup=true)
    write_amn(tmpfile, A; binary=true)
    A2 = read_amn(tmpfile)
    @test A ≈ A2
end
