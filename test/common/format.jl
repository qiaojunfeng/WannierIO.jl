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
    @test fortran_format(binary = true) isa FortranBinary
    @test fortran_format(binary = true, stream = true) isa FortranBinaryStream

    @test w90input_format() isa W90InputText
    @test w90input_format(toml = true) isa W90InputToml

    @test format_name(FortranText()) == "fortran-text"
    @test format_name(FortranBinary()) == "fortran-binary"
    @test format_name(FortranBinaryStream()) == "fortran-binary-stream"
    @test format_name(W90InputText()) == "wannier90-text"
    @test format_name(W90InputToml()) == "wannier90-toml"

    mktemp() do path, io
        write(io, "1 1 0.0\n")
        close(io)
        @test detect_fortran_format(path) isa FortranText
        @test detect_fortran_format(path; stream = true) isa FortranText
        @test detect_w90input_format(path) isa W90InputText
    end

    mktemp() do path, io
        write(io, UInt8[0x01, 0x02, 0x03, 0x04])
        close(io)
        @test detect_fortran_format(path) isa FortranBinary
        @test detect_fortran_format(path; stream = true) isa FortranBinaryStream
    end

    mktemp() do path, io
        write(io, "num_wann = 4\n")
        close(io)
        @test detect_w90input_format(path) isa W90InputToml
    end
end
