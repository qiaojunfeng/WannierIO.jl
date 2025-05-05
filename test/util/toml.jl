@testitem "write_toml" begin
    using OrderedCollections: OrderedDict

    params = (; b=2, a=1)
    # Since TOML.read returns Dict which is unordered, we need to compare literally.
    ref = """
    b = 2
    a = 1
    """
    io = IOBuffer()

    # test write_toml with NamedTuple
    WannierIO.write_toml(io, params)
    lines = String(take!(io))
    @test lines == ref

    # test write_toml with kwargs
    WannierIO.write_toml(io; params...)
    lines = String(take!(io))
    @test lines == ref

    # test write_toml with OrderedDict
    params = OrderedDict(pairs(params))
    WannierIO.write_toml(io, params)
    lines = String(take!(io))
    @test lines == ref
end
