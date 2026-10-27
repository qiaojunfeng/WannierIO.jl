@testitem "CscPack round-trip" begin
    using LazyArtifacts
    using SparseArrays: sparse

    tbdat = WannierIO.read_w90_tb_dat(artifact"Si2_valence/outputs/WS/Si2_valence_tb.dat")
    lossless = WannierIO.SparseOption(; atol = 0.0, value_type = ComplexF64, index_type = Int32)

    # lossless round-trip (ComplexF64 keeps full precision)
    csc = sparsify(tbdat.H, lossless)
    dense = densify(csc)

    @test csc isa WannierIO.CscPack{ComplexF64, Int32}
    @test size(csc) == size(tbdat.H)
    @test csc.dims == collect(size(tbdat.H))
    @test dense == tbdat.H
    @test csc[:, :, 1] == sparse(tbdat.H[:, :, 1])

    # a component axis between the matrix axes and the R axis
    r = stack((tbdat.rx, tbdat.ry, tbdat.rz); dims = 3)
    csc = sparsify(r, lossless)
    @test size(csc) == size(r)
    @test ndims(csc) == 4
    @test length(csc) == 3 * size(r, 4)
    @test densify(csc) == r
    # packed matrices follow the column-major order of the trailing axes
    iR = 5
    @test csc[:, :, 2 + 3 * (iR - 1)] == sparse(r[:, :, 2, iR])

    @test_throws ErrorException sparsify(zeros(ComplexF64, 2, 3, 4))
    @test_throws ErrorException sparsify(zeros(ComplexF64, 2, 2))
end

@testitem "SparseOperatorPack round-trip" begin
    using LazyArtifacts
    using OrderedCollections: OrderedDict

    tbdat = WannierIO.read_w90_tb_dat(artifact"Si2_valence/outputs/WS/Si2_valence_tb.dat")
    operators = OrderedDict(
        "hamiltonian" => tbdat.H,
        "berry_connection" => stack((tbdat.rx, tbdat.ry, tbdat.rz); dims = 3),
    )
    dpack = OperatorPack(tbdat.header, tbdat.lattice, tbdat.Rvectors, operators)

    # lossless round-trip (ComplexF64 keeps full precision)
    spack = sparsify(
        dpack, WannierIO.SparseOption(; atol = 0.0, value_type = ComplexF64, index_type = Int32)
    )
    dpack2 = densify(spack)

    @test spack isa WannierIO.SparseOperatorPack{Float64, Int32}
    @test spack.header == dpack.header
    @test spack.lattice == dpack.lattice
    @test n_Rvectors(spack) == n_Rvectors(dpack)
    @test n_wannier(spack) == n_wannier(dpack)
    @test spack.Rvectors == Matrix{Int32}(reduce(hcat, dpack.Rvectors))
    @test spack.operators["hamiltonian"] isa WannierIO.CscPack{ComplexF64, Int32}
    @test size(spack.operators["berry_connection"]) == size(operators["berry_connection"])
    @test dpack2.header == dpack.header
    @test dpack2.lattice == dpack.lattice
    @test dpack2.Rvectors == dpack.Rvectors
    @test dpack2.operators == dpack.operators
end

@testitem "SparseOperatorPack reduced precision" begin
    using LazyArtifacts
    using OrderedCollections: OrderedDict

    tbdat = WannierIO.read_w90_tb_dat(artifact"Si2_valence/outputs/MDRS/Si2_valence_tb.dat")
    operators = OrderedDict(
        "hamiltonian" => tbdat.H,
        "berry_connection" => stack((tbdat.rx, tbdat.ry, tbdat.rz); dims = 3),
    )
    dpack = OperatorPack(tbdat.header, tbdat.lattice, tbdat.Rvectors, operators)

    spack = sparsify(
        dpack, WannierIO.SparseOption(; atol = 0.0, value_type = Float32, index_type = Int16)
    )
    dpack32 = densify(spack)
    dpack64 = densify(spack; value_type = ComplexF64)

    @test spack isa WannierIO.SparseOperatorPack{Float32, Int16}
    @test eltype(spack.operators["hamiltonian"].nzval) == ComplexF32
    @test n_Rvectors(spack) == Int16(n_Rvectors(dpack))
    @test n_wannier(spack) == Int16(n_wannier(dpack))
    @test dpack32 isa WannierIO.OperatorPack{Float32, Int16}
    @test dpack64 isa WannierIO.OperatorPack{Float64, Int16}
    @test dpack32.lattice ≈ dpack.lattice
    @test dpack32.Rvectors == WannierIO.Vec3{Int16}.(dpack.Rvectors)
    for name in keys(operators)
        @test dpack32.operators[name] ≈ dpack.operators[name]
        @test dpack64.operators[name] ≈ ComplexF64.(dpack32.operators[name])
    end
end

@testitem "sparsify/densify atol" begin
    using LazyArtifacts
    using OrderedCollections: OrderedDict

    tbdat = WannierIO.read_w90_tb_dat(artifact"Si2_valence/outputs/MDRS/Si2_valence_tb.dat")
    operators = OrderedDict(
        "hamiltonian" => tbdat.H,
        "berry_connection" => stack((tbdat.rx, tbdat.ry, tbdat.rz); dims = 3),
    )
    dpack = OperatorPack(tbdat.header, tbdat.lattice, tbdat.Rvectors, operators)

    opt = WannierIO.SparseOption()
    spack = sparsify(dpack)

    for name in keys(dpack.operators)
        sop = spack.operators[name]
        dop = reshape(dpack.operators[name], n_wannier(dpack), n_wannier(dpack), :)
        @test all(isapprox(sop[:, :, l], dop[:, :, l]; atol = (opt.atol * 10)) for l in 1:length(sop))
        @test !all(isapprox(sop[:, :, l], dop[:, :, l]; atol = (opt.atol / 10)) for l in 1:length(sop))
    end
end
