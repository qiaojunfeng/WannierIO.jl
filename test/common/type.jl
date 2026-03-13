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
    using WannierIO: StringVec3, stringvec3, vec3

    s = "Si"
    v = [1, 2, 3]
    sv = stringvec3(s, v)

    @test (sv isa StringVec3) && (sv.first == s) && (sv.second == vec3(v))
    @test stringvec3(s, v) == sv
    @test stringvec3(s => v) == sv
    @test stringvec3(Dict(s => v)) == sv
end

@testitem "file formats" begin
    using WannierIO:
        AbstractFileFormat,
        FortranBinary,
        FortranBinaryStream,
        AbstractFortranFormat,
        FortranText,
        AbstractW90InputFormat,
        W90InputText,
        W90InputToml,
        detect_fortran_format,
        detect_w90input_format,
        fortran_format,
        format_name,
        w90input_format

    @test FortranText() isa AbstractFileFormat
    @test FortranText() isa AbstractFortranFormat
    @test W90InputToml() isa AbstractFileFormat
    @test W90InputToml() isa AbstractW90InputFormat

    @test fortran_format() isa FortranText
    @test fortran_format(binary=true) isa FortranBinary
    @test fortran_format(binary=true, stream=true) isa FortranBinaryStream

    @test w90input_format() isa W90InputText
    @test w90input_format(toml=true) isa W90InputToml

    @test format_name(FortranText()) == "fortran-text"
    @test format_name(FortranBinary()) == "fortran-binary"
    @test format_name(FortranBinaryStream()) == "fortran-binary-stream"
    @test format_name(W90InputText()) == "wannier90-text"
    @test format_name(W90InputToml()) == "wannier90-toml"

    mktemp() do path, io
        write(io, "1 1 0.0\n")
        close(io)
        @test detect_fortran_format(path) isa FortranText
        @test detect_fortran_format(path; stream=true) isa FortranText
        @test detect_w90input_format(path) isa W90InputText
    end

    mktemp() do path, io
        write(io, UInt8[0x01, 0x02, 0x03, 0x04])
        close(io)
        @test detect_fortran_format(path) isa FortranBinary
        @test detect_fortran_format(path; stream=true) isa FortranBinaryStream
    end

    mktemp() do path, io
        write(io, "num_wann = 4\n")
        close(io)
        @test detect_w90input_format(path) isa W90InputToml
    end
end
