
@testitem "isbinary file" begin
    using LazyArtifacts
    @test WannierIO.isbinary(artifact"Si2_valence/reference/binary/UNK00001.1")
    @test !WannierIO.isbinary(artifact"Si2_valence/UNK/UNK00001.1")
end

@testitem "parse_float" begin
    @test WannierIO.parse_float("1.0D-10") ≈ 1e-10
end

@testitem "parse_bool" begin
    @test WannierIO.parse_bool(".true.") == true
    @test WannierIO.parse_bool("T") == true
    @test WannierIO.parse_bool("1") == true
    @test WannierIO.parse_bool(".false") == false
end
