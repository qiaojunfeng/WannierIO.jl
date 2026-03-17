@testitem "unified tb api" begin
    using LazyArtifacts
    using HDF5
    using JLD2
    using Zarr

    tbdat = read_w90_tb_dat(artifact"Si2_valence/outputs/WS/Si2_valence_tb.dat")
    dpack0 = pack(tbdat)
    # The lattice is symmetric, let's make it asymmetric to test
    # row/column-major correctness.
    dpack0 = OperatorPack(
        dpack0.header,
        WannierIO.mat3([1.0 0.1 0.0; 0.0 1.0 0.2; 0.3 0.0 1.0]),
        dpack0.Rvectors,
        dpack0.operators,
    )

    targets = [
        (WannierIO.W90Dat(), tempname() * "_tb.dat"),
        (WannierIO.HDF5Format(), tempname() * ".h5"),
        (WannierIO.JLD2Format(), tempname() * ".jld2"),
        (WannierIO.ZarrFormat(), tempname() * ".zarr"),
        (WannierIO.ZarrZipFormat(), tempname() * ".zarr.zip"),
    ]

    for (fmt, dst) in targets
        write_w90_tb(dst, dpack0, fmt; atol=0.0, value_type=ComplexF64)
        pack2 = read_w90_tb(dst, fmt)
        pack3 = read_w90_tb(dst)

        for pack in (pack2, pack3)
            @test pack isa WannierIO.OperatorPack
            @test pack.n_wann == dpack0.n_wann
            @test pack.Rvectors == dpack0.Rvectors
            @test pack.lattice == dpack0.lattice
            for name in keys(dpack0.operators)
                @test pack.operators[name] ≈ dpack0.operators[name]
            end
        end
    end

    reduced_precision_targets = [
        (WannierIO.HDF5Format(), tempname() * ".h5"),
        (WannierIO.JLD2Format(), tempname() * ".jld2"),
        (WannierIO.ZarrFormat(), tempname() * ".zarr"),
        (WannierIO.ZarrZipFormat(), tempname() * ".zarr.zip"),
    ]

    for (fmt, dst) in reduced_precision_targets
        write_w90_tb(dst, dpack0, fmt)
        pack2 = read_w90_tb(dst, fmt)

        @test pack2 isa WannierIO.OperatorPack{Float32,Int16}
        @test pack2.n_wann == dpack0.n_wann
        @test pack2.Rvectors == dpack0.Rvectors
        @test pack2.lattice ≈ dpack0.lattice
        for name in keys(dpack0.operators)
            op2 = pack2.operators[name]
            op0 = dpack0.operators[name]
            @test isapprox(op2, op0; atol=2e-6)
            @test !isapprox(op2, op0; atol=1e-7)
        end
    end
end
