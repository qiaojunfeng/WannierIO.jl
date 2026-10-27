@testitem "OperatorPack constructor" begin
    using OrderedCollections: OrderedDict
    using WannierIO: OperatorPack, Mat3, Vec3

    header = "operator header"
    lattice = [1.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 3.0]
    Rvectors = [Vec3(0, 0, 0), Vec3(1, 0, 0)]
    operators = OrderedDict(
        "hamiltonian" => stack([ComplexF64[1 0; 0 2], ComplexF64[0 1; 1 0]]),
        "X" => stack([Float64[1 2; 3 4], Float64[5 6; 7 8]]),
        "position" => rand(ComplexF64, 2, 2, 3, 2),
    )

    pack = OperatorPack(header, lattice, Rvectors, operators)

    @test pack.header == "operator header"
    @test pack.lattice isa Mat3{Float64}
    @test pack.lattice == Mat3{Float64}(lattice)
    @test pack.Rvectors isa Vector{Vec3{Int}}
    @test pack.Rvectors == Rvectors
    @test pack.operators isa OrderedDict{String}
    @test collect(keys(pack.operators)) == ["hamiltonian", "X", "position"]
    @test eltype(pack.operators["hamiltonian"]) == ComplexF64
    @test eltype(pack.operators["X"]) == Float64
    # component axes between the matrix axes and the R axis are kept
    @test size(pack.operators["position"]) == (2, 2, 3, 2)
    @test pack.operators["position"] == operators["position"]
    @test n_Rvectors(pack) == 2
    @test n_wannier(pack) == 2
end

@testitem "OperatorPack validation" begin
    using OrderedCollections: OrderedDict
    using WannierIO: OperatorPack, Vec3

    lattice = Float64[1 0 0; 0 1 0; 0 0 1]
    Rvectors = [Vec3(0, 0, 0), Vec3(1, 0, 0)]
    build(op) = OperatorPack("bad", lattice, Rvectors, OrderedDict("O" => op))

    # last axis must be the R axis
    @test_throws ErrorException build(zeros(ComplexF64, 2, 2, 1))
    @test_throws ErrorException build(zeros(ComplexF64, 2, 2, 2, 3))
    # matrix axes must be square
    @test_throws ErrorException build(zeros(ComplexF64, 2, 3, 2))
    # a matrix without an R axis is not an operator
    @test_throws ErrorException build(zeros(ComplexF64, 2, 2))
    # all operators share n_wann
    @test_throws ErrorException OperatorPack(
        "bad", lattice, Rvectors,
        OrderedDict("A" => zeros(2, 2, 2), "B" => zeros(3, 3, 2)),
    )
end

@testitem "operator file formats" begin
    using LazyArtifacts
    using OrderedCollections: OrderedDict
    using HDF5
    using JLD2
    using Zarr

    # Raw tb.dat arrays as storage test data: the pack does no physics, so
    # the R degeneracies are irrelevant here.
    tbdat = read_w90_tb_dat(artifact"Si2_valence/outputs/WS/Si2_valence_tb.dat")
    operators = OrderedDict(
        "hamiltonian" => tbdat.H,
        "berry_connection" => stack((tbdat.rx, tbdat.ry, tbdat.rz); dims = 3),
    )
    # The lattice is symmetric, make it asymmetric to test row/column-major
    # correctness.
    lattice = [1.0 0.1 0.0; 0.0 1.0 0.2; 0.3 0.0 1.0]
    pack0 = OperatorPack(tbdat.header, lattice, tbdat.Rvectors, operators)

    targets() = [
        (WannierIO.HDF5Format(), tempname() * ".h5"),
        (WannierIO.JLD2Format(), tempname() * ".jld2"),
        (WannierIO.ZarrFormat(), tempname() * ".zarr"),
        (WannierIO.ZarrZipFormat(), tempname() * ".zarr.zip"),
    ]

    for (fmt, dst) in targets()
        write_operator(dst, pack0, fmt; atol = 0.0, value_type = ComplexF64)
        pack2 = read_operator(dst, fmt)
        pack3 = read_operator(dst)

        for pack in (pack2, pack3)
            @test pack isa WannierIO.OperatorPack
            @test pack.header == pack0.header
            @test n_wannier(pack) == n_wannier(pack0)
            @test pack.Rvectors == pack0.Rvectors
            @test pack.lattice == pack0.lattice
            @test collect(keys(pack.operators)) == ["hamiltonian", "berry_connection"]
            for name in keys(pack0.operators)
                @test size(pack.operators[name]) == size(pack0.operators[name])
                @test pack.operators[name] == pack0.operators[name]
            end
        end
    end

    for (fmt, dst) in targets()
        write_operator(dst, pack0, fmt)
        pack2 = read_operator(dst, fmt)

        @test pack2 isa WannierIO.OperatorPack{Float32, Int32}
        @test n_wannier(pack2) == n_wannier(pack0)
        @test pack2.Rvectors == pack0.Rvectors
        @test pack2.lattice ≈ pack0.lattice
        for name in keys(pack0.operators)
            op2 = pack2.operators[name]
            op0 = pack0.operators[name]
            @test size(op2) == size(op0)
            # the default `SparseOption` drops entries (and imaginary parts)
            # below `atol = 1e-6` and stores Float32: lossy per element
            @test maximum(abs, op2 - op0) < 2.0e-6
            @test op2 != op0
        end
    end
end
