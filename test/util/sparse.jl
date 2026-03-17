@testitem "CscPack round-trip" begin
    using LazyArtifacts
    using SparseArrays: sparse

    tbdat = WannierIO.read_w90_tb_dat(artifact"Si2_valence/outputs/WS/Si2_valence_tb.dat")

    # lossless round-trip (ComplexF64 keeps full precision)
    csc = sparsify(
        tbdat.H, WannierIO.SparseOption(; atol=0.0, value_type=ComplexF64, index_type=Int32)
    )
    dense = densify(csc)

    @test csc isa WannierIO.CscPack{ComplexF64,Int32}
    @test length(csc) == length(tbdat.H)
    @test WannierIO.matrix_order(csc) == size(first(tbdat.H), 1)
    @test dense == tbdat.H
    @test csc[1] == sparse(tbdat.H[1])
end

@testitem "SparseOperatorPack WS" begin
    using LazyArtifacts

    tbdat = WannierIO.read_w90_tb_dat(artifact"Si2_valence/outputs/WS/Si2_valence_tb.dat")
    dpack = pack(tbdat)

    # lossless round-trip (ComplexF64 keeps full precision)
    spack = sparsify(
        dpack, WannierIO.SparseOption(; atol=0.0, value_type=ComplexF64, index_type=Int32)
    )
    dpack2 = densify(spack)

    @test spack isa WannierIO.SparseOperatorPack{Float64,Int32}
    @test spack.header == dpack.header
    @test spack.lattice == dpack.lattice
    @test spack.n_Rvecs == dpack.n_Rvecs
    @test spack.n_wann == dpack.n_wann
    @test spack.Rvectors == Matrix{Int32}(reduce(hcat, dpack.Rvectors))
    @test spack.operators["H"] isa WannierIO.CscPack{ComplexF64,Int32}
    @test dpack2.header == dpack.header
    @test dpack2.lattice == dpack.lattice
    @test dpack2.Rvectors == dpack.Rvectors
    @test dpack2.operators == dpack.operators
end

@testitem "SparseOperatorPack MDRS" begin
    using LazyArtifacts

    tbdat = WannierIO.read_w90_tb_dat(artifact"Si2_valence/outputs/MDRS/Si2_valence_tb.dat")
    wsvec = WannierIO.read_w90_wsvec(
        artifact"Si2_valence/outputs/MDRS/Si2_valence_wsvec.dat"
    )
    dpack = pack(tbdat, wsvec)

    # reduced precision round-trip
    spack = sparsify(
        dpack, WannierIO.SparseOption(; atol=0.0, value_type=Float32, index_type=Int16)
    )
    dpack32 = densify(spack)
    dpack64 = densify(spack; value_type=ComplexF64)

    @test spack isa WannierIO.SparseOperatorPack{Float32,Int16}
    @test eltype(spack.operators["H"].nzval) == ComplexF32
    @test spack.n_Rvecs == Int16(dpack.n_Rvecs)
    @test spack.n_wann == Int16(dpack.n_wann)
    @test dpack32 isa WannierIO.OperatorPack{Float32,Int16}
    @test dpack64 isa WannierIO.OperatorPack{Float64,Int16}
    @test dpack32.lattice ≈ dpack.lattice
    @test dpack32.Rvectors == WannierIO.Vec3{Int16}.(dpack.Rvectors)
    @test dpack32.operators["H"] ≈ dpack.operators["H"]
    @test dpack64.operators["H"] ≈ [ComplexF64.(O) for O in dpack32.operators["H"]]
end

@testitem "sparsify/densify atol" begin
    using LazyArtifacts

    tbdat = WannierIO.read_w90_tb_dat(artifact"Si2_valence/outputs/MDRS/Si2_valence_tb.dat")
    dpack = pack(tbdat)

    # reduced precision round-trip
    opt = WannierIO.SparseOption()
    spack = sparsify(dpack)

    for name in keys(dpack.operators)
        sop = spack.operators[name]
        dop = dpack.operators[name]
        @test all(isapprox.(sop, dop; atol=(opt.atol * 10)))
        @test !all(isapprox.(sop, dop; atol=(opt.atol / 10)))
    end
end
