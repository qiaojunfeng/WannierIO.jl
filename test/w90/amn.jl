@testset "read/write amn" begin
    A = read_amn(artifact"Si2_valence/si2.amn")

    tmpfile = tempname(; cleanup=true)
    write_amn(tmpfile, A)
    A2 = read_amn(tmpfile)
    @test A ≈ A2
end

@testset "read/write amn binary" begin
    A = read_amn(artifact"Si2_valence/si2.amn")
    A1 = read_amn(artifact"Si2_valence/reference/si2.binary.amn")
    @test A ≈ A1

    tmpfile = tempname(; cleanup=true)
    write_amn(tmpfile, A; binary=true)
    A2 = read_amn(tmpfile)
    @test A ≈ A2
end
