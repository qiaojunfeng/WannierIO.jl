@testitem "unified tb api" begin
    using HDF5
    using JLD2
    using LazyArtifacts
    using Zarr

    tbdat = read_w90_tb_dat(artifact"Si2_valence/outputs/WS/Si2_valence_tb.dat")
    dpack0 = pack(tbdat)

    targets = [
        (WannierIO.W90Dat(), tempname() * "_tb.dat"),
        (WannierIO.HDF5Format(), tempname() * ".h5"),
        (WannierIO.JLD2Format(), tempname() * ".jld2"),
        (WannierIO.ZarrFormat(), tempname() * ".zarr"),
    ]

    for (fmt, dst) in targets
        write_w90_tb(dst, dpack0, fmt; atol=0.0, value_type=ComplexF64)
        pack2 = read_w90_tb(dst, fmt)
        pack3 = read_w90_tb(dst)

        for pack in (pack2, pack3)
            @test pack isa WannierIO.OperatorPack
            @test pack.n_wann == dpack0.n_wann
            @test pack.Rvectors == dpack0.Rvectors
            for name in keys(dpack0.operators)
                @test pack.operators[name] ≈ dpack0.operators[name]
            end
        end
    end

    reduced_precision_targets = [
        (WannierIO.HDF5Format(), tempname() * ".h5"),
        (WannierIO.JLD2Format(), tempname() * ".jld2"),
        (WannierIO.ZarrFormat(), tempname() * ".zarr"),
    ]

    for (fmt, dst) in reduced_precision_targets
        write_w90_tb(dst, dpack0, fmt; atol=0.0, value_type=Float32, index_type=Int16)
        pack2 = read_w90_tb(dst, fmt)
        @test pack2 isa WannierIO.OperatorPack{Float32,Int16}
        @test pack2.lattice ≈ dpack0.lattice
        @test pack2.operators["H"] ≈ dpack0.operators["H"]
    end
end
