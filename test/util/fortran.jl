
@testset "isbinary file" begin
    @test WannierIO.isbinary(artifacts"Si2/UNK/binary/UNK00001.1")
    @test !WannierIO.isbinary(artifacts"Si2/UNK/textual/UNK00001.1")
end

@testset "parse_float" begin
    @test WannierIO.parse_float("1.0D-10") ≈ 1e-10
end

@testset "parse_bool" begin
    @test WannierIO.parse_bool(".true.") == true
    @test WannierIO.parse_bool("T") == true
    @test WannierIO.parse_bool("1") == true
    @test WannierIO.parse_bool(".false") == false
end