@testitem "CscPack round-trip" begin
    using LazyArtifacts
    using SparseArrays: sparse

    tbdat = WannierIO.read_w90_tb_dat(artifact"Si2_valence/outputs/WS/Si2_valence_tb.dat")

    # lossless by default, at the element type of the operator
    csc = CscPack(tbdat.H)
    @test csc isa CscPack{ComplexF64, Int32}
    @test size(csc) == size(tbdat.H)
    @test csc.dims == collect(size(tbdat.H))
    @test Array(csc) == tbdat.H
    @test Array{ComplexF32}(csc) ≈ tbdat.H
    @test csc[:, :, 1] == sparse(tbdat.H[:, :, 1])

    # a component axis between the matrix axes and the R axis
    r = stack((tbdat.rx, tbdat.ry, tbdat.rz); dims = 3)
    csc = CscPack(r)
    @test size(csc) == size(r)
    @test ndims(csc) == 4
    @test length(csc) == 3 * size(r, 4)
    @test Array(csc) == r
    # packed matrices follow the column-major order of the trailing axes
    iR = 5
    @test csc[:, :, 2 + 3 * (iR - 1)] == sparse(r[:, :, 2, iR])

    @test CscPack{ComplexF32, Int16}(r) isa CscPack{ComplexF32, Int16}
    @test CscPack{ComplexF32}(r) isa CscPack{ComplexF32, Int32}

    @test_throws ArgumentError CscPack(zeros(ComplexF64, 2, 2, 0))
    @test_throws ErrorException CscPack(zeros(ComplexF64, 2, 3, 4))
    @test_throws ErrorException CscPack(zeros(ComplexF64, 2, 2))
    @test_throws ArgumentError CscPack(r; atol = -1)
end

@testitem "CscPack atol" begin
    op = zeros(ComplexF64, 3, 3, 2)
    op[1, 1, 1] = 1.0e-300
    op[2, 3, 1] = 2.0 + 1.0e-9im
    op[3, 2, 2] = 0.0 - 1.0e-300im

    # atol = 0 drops exactly the zeros
    csc = CscPack(op)
    @test csc.nzptr[end] - 1 == 3
    @test Array(csc) == op

    # an entry is dropped when both parts are at most atol; a kept entry
    # loses its part that is at most atol
    csc = CscPack(op; atol = 1.0e-9)
    @test csc.nzptr[end] - 1 == 1
    @test Array(csc)[2, 3, 1] === 2.0 + 0.0im

    # a real value type takes an operator whose remaining imaginary parts vanish
    @test Array(CscPack{Float64}(op; atol = 1.0e-9))[2, 3, 1] === 2.0
    @test_throws InexactError CscPack{Float64}(op)
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

    # lossless by default
    spack = SparseOperatorPack(dpack)
    dpack2 = OperatorPack(spack)

    @test spack isa SparseOperatorPack{Float64, Int32}
    @test spack.header == dpack.header
    @test spack.lattice == dpack.lattice
    @test n_Rvectors(spack) == n_Rvectors(dpack)
    @test n_wannier(spack) == n_wannier(dpack)
    @test spack.Rvectors == Matrix{Int32}(reduce(hcat, dpack.Rvectors))
    @test spack.operators["hamiltonian"] isa CscPack{ComplexF64, Int32}
    @test size(spack.operators["berry_connection"]) == size(operators["berry_connection"])
    @test dpack2 isa OperatorPack{Float64, Int32}
    @test dpack2.header == dpack.header
    @test dpack2.lattice == dpack.lattice
    @test dpack2.Rvectors == dpack.Rvectors
    @test dpack2.operators == dpack.operators

    @test_throws ArgumentError SparseOperatorPack(
        OperatorPack("empty", tbdat.lattice, tbdat.Rvectors, OrderedDict{String, Array{Float64, 3}}())
    )
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

    spack = SparseOperatorPack{Float32, Int16}(dpack)
    dpack32 = OperatorPack(spack)
    dpack64 = OperatorPack{Float64}(spack)

    @test SparseOperatorPack{Float32}(dpack) isa SparseOperatorPack{Float32, Int32}
    @test spack isa SparseOperatorPack{Float32, Int16}
    @test eltype(spack.operators["hamiltonian"].nzval) == ComplexF32
    @test n_Rvectors(spack) == n_Rvectors(dpack)
    @test n_wannier(spack) == n_wannier(dpack)
    @test dpack32 isa OperatorPack{Float32, Int16}
    @test dpack64 isa OperatorPack{Float64, Int16}
    @test dpack32.lattice == dpack.lattice
    @test spack.lattice isa WannierIO.Mat3{Float64}
    @test dpack32.Rvectors == WannierIO.Vec3{Int16}.(dpack.Rvectors)
    for name in keys(operators)
        # Float32 rounds each element, it drops nothing
        @test dpack32.operators[name] ≈ dpack.operators[name]
        @test dpack64.operators[name] == ComplexF64.(dpack32.operators[name])
    end
end

@testitem "SparseOperatorPack atol" begin
    using LazyArtifacts
    using OrderedCollections: OrderedDict

    tbdat = WannierIO.read_w90_tb_dat(artifact"Si2_valence/outputs/MDRS/Si2_valence_tb.dat")
    operators = OrderedDict(
        "hamiltonian" => tbdat.H,
        "berry_connection" => stack((tbdat.rx, tbdat.ry, tbdat.rz); dims = 3),
    )
    dpack = OperatorPack(tbdat.header, tbdat.lattice, tbdat.Rvectors, operators)

    atol = 1.0e-6
    lossless = SparseOperatorPack(dpack)
    spack = SparseOperatorPack(dpack; atol)
    back = OperatorPack(spack)
    for name in keys(dpack.operators)
        # each part changes by at most atol
        @test maximum(abs, back.operators[name] - dpack.operators[name]) <= sqrt(2) * atol
    end
    n_stored(p) = sum(csc -> length(csc.nzval), values(p.operators))
    @test n_stored(spack) < n_stored(lossless)
end

@testitem "SparseOperatorPack of a real-valued operator" begin
    using OrderedCollections: OrderedDict
    using WannierIO: Vec3

    lattice = Float64[1 0 0; 0 1 0; 0 0 1]
    Rvectors = [Vec3(0, 0, 0), Vec3(1, 0, 0)]
    H = ComplexF64.(stack([Float64[1 0; 0 2], Float64[0 3; 3 0]]))
    r = rand(ComplexF64, 2, 2, 3, 2)
    dpack = OperatorPack("real", lattice, Rvectors, OrderedDict("H" => H, "r" => r))

    # an operator with vanishing imaginary parts is stored, and read back, real
    spack = SparseOperatorPack(dpack)
    @test spack.operators["H"] isa CscPack{Float64}
    @test spack.operators["r"] isa CscPack{ComplexF64}
    back = OperatorPack(spack)
    @test back.operators["H"] isa Array{Float64}
    @test back.operators["H"] == H
    @test back.operators["r"] == r
end
