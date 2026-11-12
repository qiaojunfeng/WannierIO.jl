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

    # the operators set the value precision, the lattice stays Float64
    lattice_frac = [3.77943445 -1.88971723 0.0; 0.0 3.27308625 0.0; 0.0 0.0 6.05751038]
    pack32 = OperatorPack(
        header, lattice_frac, Rvectors, OrderedDict("H" => rand(ComplexF32, 2, 2, 2))
    )
    @test pack32 isa OperatorPack{Float32}
    @test pack32.lattice isa Mat3{Float64}
    @test pack32.lattice == lattice_frac
    # mixed precisions are promoted
    mixed = OrderedDict("A" => rand(Float32, 2, 2, 2), "B" => rand(ComplexF64, 2, 2, 2))
    @test OperatorPack(header, lattice, Rvectors, mixed) isa OperatorPack{Float64}
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

    # a file holds a SparseOperatorPack and reads back as written
    for spack0 in (SparseOperatorPack(pack0), SparseOperatorPack{Float32}(pack0; atol = 1.0e-6))
        for (fmt, dst) in targets()
            write_operator(dst, spack0, fmt)
            for spack in (read_operator(dst, fmt), read_operator(dst))
                @test spack isa typeof(spack0)
                @test spack.header == spack0.header
                @test spack.Rvectors == spack0.Rvectors
                # the lattice is kept in Float64 whatever the precision of the operators
                @test spack.lattice == pack0.lattice
                @test collect(keys(spack.operators)) == ["hamiltonian", "berry_connection"]
                for (name, csc0) in spack0.operators
                    csc = spack.operators[name]
                    @test typeof(csc) == typeof(csc0)
                    for field in fieldnames(CscPack)
                        @test getfield(csc, field) == getfield(csc0, field)
                    end
                end
            end
        end
    end

    # lossless storage gives back the dense pack
    for (fmt, dst) in targets()
        write_operator(dst, SparseOperatorPack(pack0), fmt)
        pack = OperatorPack(read_operator(dst))
        @test pack.operators == pack0.operators
    end

    # a dense pack must be converted first
    for (fmt, dst) in targets()
        @test_throws ArgumentError write_operator(dst, pack0, fmt)
        @test_throws ArgumentError write_operator(dst, pack0)
    end
end
