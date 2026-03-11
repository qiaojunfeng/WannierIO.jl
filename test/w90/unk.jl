
@testitem "read/write unk" begin
    using LazyArtifacts
    unk = read_unk(artifact"Si2_valence/UNK/UNK00001.1")
    ik = unk.ik
    ψ = unk.Ψ
    @test ik == 1
    @test size(ψ) == (13, 13, 13, 4, 1)
    ref_ψ234 = ComplexF64[
        -0.82739368086-1.2362766865im; -1.1707908501+0.11631095198im; -0.25943725378+2.3452525646im; 1.090429898-1.8525861633im;;
    ]
    @test ψ[2, 3, 4, :, :] == ref_ψ234

    tmpfile = tempname(; cleanup=true)
    write_unk(tmpfile, ik, ψ)
    unk2 = read_unk(tmpfile)
    ik2 = unk2.ik
    ψ2 = unk2.Ψ

    @test ik ≈ ik2
    @test ψ ≈ ψ2
end

@testitem "read/write unk binary" begin
    using LazyArtifacts
    unk = read_unk(artifact"Si2_valence/UNK/UNK00001.1")
    ik = unk.ik
    ψ = unk.Ψ
    unk1 = read_unk(artifact"Si2_valence/outputs/binary/UNK00001.1")
    ik1 = unk1.ik
    ψ1 = unk1.Ψ
    @test ik == ik1
    @test ψ ≈ ψ1

    tmpfile = tempname(; cleanup=true)
    write_unk(tmpfile, ik, ψ; binary=true)
    unk2 = read_unk(tmpfile)
    ik2 = unk2.ik
    ψ2 = unk2.Ψ

    @test ik == ik2
    @test ψ ≈ ψ2
end
