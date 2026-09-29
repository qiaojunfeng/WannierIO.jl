@testitem "read nnkp" begin
    using LazyArtifacts
    nnkp = read_nnkp(artifact"Si2_valence/outputs/Si2_valence.nnkp")

    WRITE_TOML = false
    WRITE_TOML && write_nnkp("/tmp/Si2_valence.nnkp.toml", nnkp, WannierIO.W90InputToml())

    test_data = read_nnkp(artifact"Si2_valence/outputs/Si2_valence.nnkp.toml")
    @test nnkp isa Nnkp
    @test nnkp == test_data
    @test length(nnkp.projections) == 4
    @test isnothing(nnkp.spinor_projections)
    @test isnothing(nnkp.auto_projections)
    @test isempty(nnkp.exclude_bands)
end

@testitem "read/write nnkp" begin
    using LazyArtifacts
    nnkp = read_nnkp(artifact"Si2_valence/outputs/Si2_valence.nnkp")
    tmpfile = tempname(; cleanup = true)
    write_nnkp(tmpfile, nnkp)

    nnkp2 = read_nnkp(tmpfile)
    @test nnkp == nnkp2
end

@testitem "read/write nnkp toml" begin
    using LazyArtifacts
    # Note that this requires https://github.com/JuliaLang/julia/pull/57584
    if VERSION > v"1.11.4"
        nnkp = read_nnkp(artifact"Si2_valence/outputs/Si2_valence.nnkp.toml")

        tmpfile = tempname(; cleanup = true)
        write_nnkp(tmpfile, nnkp, WannierIO.W90InputToml())

        nnkp2 = read_nnkp(tmpfile)
        @test nnkp == nnkp2
    end
end

@testitem "read spinor_projections" begin
    using LazyArtifacts
    nnkp = read_nnkp(artifact"Fe_soc/outputs/Fe.nnkp")

    @test isnothing(nnkp.projections)
    sp = nnkp.spinor_projections
    @test sp isa AbstractVector{<:WannierIO.SpinorHydrogenOrbital}
    @test length(sp) == 16
    # spin alternates +1 / -1 for consecutive up/down partners
    @test sp[1].spin == 1
    @test sp[2].spin == -1
    @test sp[1].spin_qaxis == [0.0, 0.0, 1.0]
    # the exclude_bands block is parsed, not skipped
    @test nnkp.exclude_bands == 1:8
end

@testitem "read/write spinor_projections" begin
    using LazyArtifacts
    nnkp = read_nnkp(artifact"Fe_soc/outputs/Fe.nnkp")
    tmpfile = tempname(; cleanup = true)
    write_nnkp(tmpfile, nnkp)

    nnkp2 = read_nnkp(tmpfile)
    @test nnkp == nnkp2
end

@testitem "read/write spinor_projections toml" begin
    using LazyArtifacts
    # Note that this requires https://github.com/JuliaLang/julia/pull/57584
    if VERSION > v"1.11.4"
        nnkp = read_nnkp(artifact"Fe_soc/outputs/Fe.nnkp")

        tmpfile = tempname(; cleanup = true)
        write_nnkp(tmpfile, nnkp, WannierIO.W90InputToml())

        nnkp2 = read_nnkp(tmpfile)
        @test nnkp == nnkp2
    end
end

@testitem "read/write empty spinor_projections (auto + spinor)" begin
    using LinearAlgebra: I
    # Wannier90 writes an empty `spinor_projections` block together with an
    # `auto_projections` block when using automatic projections in a spinor
    # (noncollinear) calculation.
    params = Nnkp(;
        lattice = WannierIO.Mat3(Matrix(1.0I, 3, 3)),
        recip_lattice = WannierIO.Mat3(Matrix(1.0I, 3, 3)),
        kpoints = [[0.0, 0.0, 0.0]],
        kpb_k = reshape([1], 1, 1),
        kpb_G = reshape([WannierIO.Vec3(0, 0, 0)], 1, 1),
        spinor_projections = WannierIO.SpinorHydrogenOrbital[],
        auto_projections = 4,
    )

    # text round-trip: an empty block stays distinct from an absent one
    tmp = tempname(; cleanup = true)
    write_nnkp(tmp, params)
    p = read_nnkp(tmp)
    @test p == params
    @test isempty(p.spinor_projections)
    @test isnothing(p.projections)
    @test p.auto_projections == 4

    # TOML round-trip (the empty list must stay typed as SpinorHydrogenOrbital)
    if VERSION > v"1.11.4"
        tmp2 = tempname(; cleanup = true)
        write_nnkp(tmp2, params, WannierIO.W90InputToml())
        p2 = read_nnkp(tmp2)
        @test p2 == params
        @test isempty(p2.spinor_projections)
        @test isnothing(p2.projections)
    end
end

@testitem "write nnkp validates projection element type" begin
    using LazyArtifacts
    # Read a real spinor nnkp to get genuine SpinorHydrogenOrbital data
    nnkp = read_nnkp(artifact"Fe_soc/outputs/Fe.nnkp")

    # Spinor orbitals as plain `projections` (or plain orbitals as
    # `spinor_projections`) must be rejected rather than written to the wrong block.
    blocks = (; nnkp.lattice, nnkp.recip_lattice, nnkp.kpoints, nnkp.kpb_k, nnkp.kpb_G)
    @test_throws ArgumentError Nnkp(; blocks..., projections = nnkp.spinor_projections)
    plain = [WannierIO.HydrogenOrbital(; (k => getfield(o, k) for k in fieldnames(WannierIO.HydrogenOrbital))...) for o in nnkp.spinor_projections]
    @test_throws ArgumentError Nnkp(; blocks..., spinor_projections = plain)
end

@testitem "read auto_projections" begin
    using LazyArtifacts
    nnkp = read_nnkp(artifact"SnSe2/outputs/SnSe2.nnkp")
    @test nnkp.auto_projections == 12
    # an empty `projections` block accompanies `auto_projections`
    @test isempty(nnkp.projections)
    @test nnkp.exclude_bands == 1:5

    tmpfile = tempname(; cleanup = true)
    write_nnkp(tmpfile, nnkp, WannierIO.W90InputToml())
    nnkp2 = read_nnkp(tmpfile)
    @test nnkp2 == nnkp
end
