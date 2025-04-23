@testitem "parse_vector" begin
    io = IOBuffer("""1  2  3  4  5  6  7  8  9  10
    11 12 13 14 15 16 17 18 19 20
    21 22 23""")
    vec = WannierIO.parse_vector(io, Int, 23)
    @test vec == collect(1:23)
    close(io)

    # different number of elements per line
    io = IOBuffer("""1  2  3  4  5
    6
    7 8 9""")
    vec = WannierIO.parse_vector(io, Int, 9)
    @test vec == collect(1:9)
    close(io)
end

@testitem "format_indices" begin
    res = WannierIO.format_indices([1, 2, 5, 8, 9, 10])
    ref = "1-2, 5, 8-10"
    @test res == ref

    res = WannierIO.format_indices(1:2)
    ref = "1-2"
    @test res == ref
end
