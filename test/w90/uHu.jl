@testitem "read/write uHu" begin
    using LazyArtifacts
    uHu = read_uHu(artifact"Fe_soc/test/Fe.uHu.fmt")

    ref_uHu_end = ComplexF64[
        22.832600673-6.6071280391e-16im 0.00020289010265+0.00033152537834im 0.00015694794962+0.00019375661933im -0.011251642778-0.09153385574im
        0.00020289010265-0.00033152537835im 26.088879829-6.3490879221e-16im -1.0614925817+0.14032143069im -0.0063702801404-0.0029107828891im
        0.00015694794962-0.00019375661933im -1.0614925817-0.14032143069im 24.575166569-1.0842021725e-19im -0.0036277162277-0.0022592847918im
        -0.011251642778+0.09153385574im -0.0063702801404+0.0029107828891im -0.0036277162277+0.0022592847918im 28.464375953-1.1301723446e-15im
    ]

    @test length(uHu) == 8
    @test size(uHu[1]) == (12, 12)
    @test size(uHu[1][1, 1]) == (4, 4)
    @test uHu[end][end, end][1:4, 1:4] ≈ ref_uHu_end

    tmpfile = tempname(; cleanup=true)
    write_uHu(tmpfile, uHu; binary=false)

    uHu1 = read_uHu(tmpfile)
    @test uHu ≈ uHu1
end

@testitem "read/write uHu binary" begin
    using LazyArtifacts
    uHu = read_uHu(artifact"Fe_soc/test/Fe.uHu")
    uHu1 = read_uHu(artifact"Fe_soc/test/Fe.uHu.fmt")
    @test uHu ≈ uHu1

    tmpfile = tempname(; cleanup=true)
    write_uHu(tmpfile, uHu; binary=true)
    uHu2 = read_uHu(tmpfile)
    @test uHu ≈ uHu2
end
