
@testset "read/write eig" begin
    E = read_eig(artifact"Si2_valence/si2.eig")
    @test length(E) == 216
    @test length(E[1]) == 4

    ref_E2 = [-5.463742712458, 3.764318343758, 5.709948162634, 5.709948162651]
    @test E[2] ≈ ref_E2

    tmpfile = tempname(; cleanup=true)
    write_eig(tmpfile, E)
    E1 = read_eig(tmpfile)

    @test E ≈ E1
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
