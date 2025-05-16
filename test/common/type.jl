@testitem "vec3 mat3" begin
    using WannierIO: Vec3, vec3, Mat3, mat3

    v0 = [1, 2, 3]
    v1 = Vec3{Int}(v0)

    @test vec3(v1) === v1
    v2 = vec3(v0)
    @test (v2 isa Vec3) && (v2 == v1)

    m0 = [1 2 3; 4 5 6; 7 8 9]
    m1 = Mat3{Int}(m0)

    @test mat3(m1) === m1
    m2 = mat3(m0)
    @test (m2 isa Mat3) && (m2 == m1)
    # To Vector{Vector{Int}}
    m3 = collect.(eachcol(m0))
    m4 = mat3(m3)
    @test (m4 isa Mat3) && (m4 == m1)

    @test vec3(m1) == Vec3{Vec3{Int}}(m3[1], m3[2], m3[3])
end

@testitem "symbolvec3" begin
    using WannierIO: SymbolVec3, symbolvec3, vec3

    s = :Si
    v = [1, 2, 3]
    sv = symbolvec3(s, v)

    @test (sv isa SymbolVec3) && (sv.first == s) && (sv.second == vec3(v))
    @test symbolvec3(string(s), v) == sv
    @test symbolvec3(s => v) == sv
    @test symbolvec3(Dict(s => v)) == sv
end
