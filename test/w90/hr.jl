@testitem "read hr ws" begin
    using LazyArtifacts
    hrdat = read_w90_hr_dat(artifact"Si2_valence/outputs/WS/Si2_valence_hr.dat")
    # just some simple tests
    @test length(hrdat.Rvectors) == 279
    @test hrdat.Rvectors[1] == [-4, 0, 2]
    @test length(hrdat.H) == 279
    H1 = ComplexF64[
        0.000805+0.0im -0.000431-0.0im -0.000137+0.0im -0.000431+0.0im
        -0.000431-0.0im 0.000805+0.0im -0.000137+0.0im -0.000431-0.0im
        -0.000299+0.0im -0.000298+0.0im 0.000534+0.0im -0.000298-0.0im
        -0.000431+0.0im -0.000431+0.0im -0.000137-0.0im 0.000805-0.0im
    ]
    @test hrdat.H[1] ≈ H1
end

@testitem "read hr mdrs" begin
    using LazyArtifacts
    # The two files are identical
    hrdat = read_w90_hr_dat(artifact"Si2_valence/outputs/WS/Si2_valence_hr.dat")
    hrdat_ws = read_w90_hr_dat(artifact"Si2_valence/outputs/WS/Si2_valence_hr.dat")

    for p in propertynames(hrdat)
        p == :header && continue
        @test getfield(hrdat_ws, p) ≈ getfield(hrdat, p)
    end
end

@testitem "write hr" begin
    using LazyArtifacts
    hrdat = read_w90_hr_dat(artifact"Si2_valence/outputs/WS/Si2_valence_hr.dat")

    tmpfile = tempname(; cleanup=true)
    write_w90_hr_dat(tmpfile, hrdat)
    hrdat2 = read_w90_hr_dat(tmpfile)

    @test propertynames(hrdat) == propertynames(hrdat2)
    for name in propertynames(hrdat)
        name == :header && continue
        @test getfield(hrdat2, name) == getfield(hrdat, name)
    end
end

@testitem "pack hr applies Rdegens" begin
    using WannierIO: Vec3

    Rvectors = [Vec3(0, 0, 0), Vec3(1, 0, 0)]
    Rdegens = [2, 4]
    H = [
        ComplexF64[2 0; 0 6],
        ComplexF64[4 0; 0 8],
    ]
    hrdat = WannierIO.HrDat("hr test", Rvectors, Rdegens, H)

    p = pack(hrdat)

    @test p.operators["H"][1] == ComplexF64[1 0; 0 3]
    @test p.operators["H"][2] == ComplexF64[1 0; 0 2]
end

@testitem "unified hr api" begin
    using LazyArtifacts
    using HDF5
    using JLD2
    using Zarr

    tbdat = read_w90_tb_dat(artifact"Si2_valence/outputs/WS/Si2_valence_tb.dat")
    dpack0 = pack(tbdat)
    hrdat0 = WannierIO.HrDat(dpack0)

    # Native W90Dat format
    w90_file = tempname() * "_hr.dat"
    write_w90_hr(w90_file, hrdat0)
    pack_w90 = read_w90_hr(w90_file)
    @test pack_w90 isa WannierIO.OperatorPack
    @test pack_w90.n_wann == dpack0.n_wann
    @test pack_w90.Rvectors == dpack0.Rvectors
    @test collect(keys(pack_w90.operators)) == ["H"]
    @test all(isnan, pack_w90.lattice)
    @test isapprox(pack_w90.operators["H"], dpack0.operators["H"]; atol=2e-5)

    # Backend formats via read_operators/write_operators
    targets = [
        (WannierIO.HDF5Format(), tempname() * ".h5"),
        (WannierIO.JLD2Format(), tempname() * ".jld2"),
        (WannierIO.ZarrFormat(), tempname() * ".zarr"),
        (WannierIO.ZarrZipFormat(), tempname() * ".zarr.zip"),
    ]

    for (fmt, dst) in targets
        write_operators(dst, dpack0, fmt; atol=0.0, value_type=ComplexF64)
        pack2 = read_operators(dst, fmt)
        
        @test pack2 isa WannierIO.OperatorPack
        @test pack2.n_wann == dpack0.n_wann
        @test pack2.Rvectors == dpack0.Rvectors
        @test collect(keys(pack2.operators)) == ["H", "r_x", "r_y", "r_z"]
        @test isapprox(pack2.operators["H"], dpack0.operators["H"]; atol=0.0)
    end
end
