
@testset "read/write eig" begin
    E = read_eig(artifact"Si2_valence/si2.eig")

    tmpfile = tempname(; cleanup=true)
    write_eig(tmpfile, E)
    E2 = read_eig(tmpfile)

    @test E ≈ E2
end

@testset "read/write eig binary" begin
    E = read_eig(artifact"Si2_valence/si2.eig")
    E1 = read_eig(artifact"Si2_valence/reference/binary/si2.eig")
    @test E ≈ E1

    tmpfile = tempname(; cleanup=true)
    write_eig(tmpfile, E; binary=true)
    E2 = read_eig(tmpfile)
    @test E ≈ E2
end
