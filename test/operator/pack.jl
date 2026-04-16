@testitem "OperatorPack constructor" begin
    using OrderedCollections: OrderedDict
    using WannierIO: OperatorPack, Mat3, Vec3

    header = "operator header"
    lattice = [1.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 3.0]
    Rvectors = [Vec3(0, 0, 0), Vec3(1, 0, 0)]
    operators = OrderedDict(
        "H" => stack([ComplexF64[1 0; 0 2], ComplexF64[0 1; 1 0]]),
        "X" => stack([Float64[1 2; 3 4], Float64[5 6; 7 8]]),
    )

    pack = OperatorPack(header, lattice, Rvectors, operators)

    @test pack.header == "operator header"
    @test pack.lattice isa Mat3{Float64}
    @test pack.lattice == Mat3{Float64}(lattice)
    @test pack.Rvectors isa Vector{Vec3{Int}}
    @test pack.Rvectors == Rvectors
    @test pack.operators isa OrderedDict{String}
    @test collect(keys(pack.operators)) == ["H", "X"]
    @test eltype(eltype(pack.operators["H"])) == ComplexF64
    @test eltype(eltype(pack.operators["X"])) == Float64
    @test n_Rvectors(pack) == 2
    @test n_wannier(pack) == 2
end

@testitem "OperatorPack validation" begin
    using OrderedCollections: OrderedDict
    using WannierIO: OperatorPack, Vec3

    lattice = Float64[1 0 0; 0 1 0; 0 0 1]
    Rvectors = [Vec3(0, 0, 0), Vec3(1, 0, 0)]

    bad_length_ops = OrderedDict("H" => [ComplexF64[1 0; 0 1]])
    @test_throws ErrorException OperatorPack(
        "bad length", lattice, Rvectors, bad_length_ops
    )

    bad_shape_ops = OrderedDict("H" => [ComplexF64[1 0 0; 0 1 0], ComplexF64[1 0 0; 0 1 0]])
    @test_throws ErrorException OperatorPack("bad shape", lattice, Rvectors, bad_shape_ops)
end
