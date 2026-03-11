@testitem "strip_comment" begin
    using WannierIO: strip_comment

    line = " 1 2 3 ! this is a comment"
    res = strip_comment(line; spaces=true)
    @test res == " 1 2 3 "

    line = " 1 2 3 # this is a comment"
    res = strip_comment(line; spaces=true)
    @test res == " 1 2 3 "

    line = " 1 2 3 "
    res = strip_comment(line; spaces=true)
    @test res == " 1 2 3 "

    line = " # 1 2 3 "
    res = strip_comment(line; spaces=true)
    @test res == " "

    line = " 1 2 3 ! this is a comment # with a hash"
    res = strip_comment(line)
    @test res == "1 2 3"
end

@testitem "parse_vector" begin
    using WannierIO: parse_vector

    io = IOBuffer("""1  2  3  4  5  6  7  8  9  10
    11 12 13 14 15 16 17 18 19 20
    21 22 23""")
    vec = parse_vector(io, Int, 23)
    @test vec == collect(1:23)
    close(io)

    # different number of elements per line
    io = IOBuffer("""1  2  3  4  5
    6
    7 8 9""")
    vec = parse_vector(io, Int, 9)
    @test vec == collect(1:9)
    close(io)

    line = "1 2 3 4 5"
    vec = parse_vector(line, Int)
    @test vec == [1, 2, 3, 4, 5]
end

@testitem "format_indices" begin
    res = WannierIO.format_indices([1, 2, 5, 8, 9, 10])
    ref = "1-2, 5, 8-10"
    @test res == ref

    res = WannierIO.format_indices(1:2)
    ref = "1-2"
    @test res == ref
end

@testitem "parse_indices" begin
    s = "1-2, 5,8 -10"
    ref = [1, 2, 5, 8, 9, 10]
    idxs = WannierIO.parse_indices(s)
    @test idxs == ref

    s = "1-2 5 8 -10"
    idxs = WannierIO.parse_indices(s)
    @test idxs == ref
end

@testitem "nextline" begin
    using WannierIO: nextline

    # Test 1: single non-empty line
    io = IOBuffer("hello")
    @test nextline(io) == "hello"
    close(io)

    # Test 2: multiple empty lines followed by non-empty line
    io = IOBuffer(join(["", "", "data"], "\n"))
    @test nextline(io) == "data"
    close(io)

    # Test 3: lines with only whitespace should be ignored
    io = IOBuffer(join(["   ", "\t", "  \t  ", "content"], "\n"))
    @test nextline(io) == "content"
    close(io)

    # Test 4: non-empty line with leading/trailing whitespace
    io = IOBuffer("   text with spaces   ")
    @test nextline(io) == "text with spaces"
    close(io)

    # Test 5: read multiple non-empty lines sequentially
    io = IOBuffer(join(["line1", "", "line2", "", "line3"], "\n"))
    @test nextline(io) == "line1"
    @test nextline(io) == "line2"
    @test nextline(io) == "line3"
    close(io)

    # Test 6: EOF without finding non-empty line should return empty string
    io = IOBuffer("")
    @test nextline(io) === ""
    close(io)

    # Test 7: only empty lines should return empty string
    io = IOBuffer(join(["", "  ", "\t"], "\n"))
    @test nextline(io) === ""
    close(io)

    # Test 8: removes comments
    io = IOBuffer("data ! this is a comment")
    @test nextline(io) == "data"
    close(io)

    # Test 9: comment with hash
    io = IOBuffer("value # comment here")
    @test nextline(io) == "value"
    close(io)

    # Test 10: comment with multiple empty lines
    io = IOBuffer(join(["", "", "data ! ignored comment", "", "", "next"], "\n"))
    @test nextline(io) == "data"
    @test nextline(io) == "next"
    close(io)

    # Test 11: keep comment
    io = IOBuffer("data ! this is a comment")
    @test nextline(io; comment=true) == "data ! this is a comment"
    close(io)

    # Test 12: lowercase option
    io = IOBuffer("HELLO WORLD")
    @test nextline(io; lower=true) == "hello world"
    close(io)

    # Test 13: lowercase with empty lines
    io = IOBuffer(join(["", "", "MIXED Case"], "\n"))
    @test nextline(io; lower=true) == "mixed case"
    close(io)

    # Test 14: both comment and lower
    io = IOBuffer("DATA VALUE ! comment")
    @test nextline(io; comment=true, lower=true) == "data value ! comment"
    close(io)
end
