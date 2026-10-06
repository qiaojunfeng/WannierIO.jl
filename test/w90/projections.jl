@testitem "parse_projections atom labels and order" begin
    lattice = [0 1 1; 1 0 1; 1 1 0] * 2.715
    atoms = ["Si" => [0.0, 0.0, 0.0], "Ga" => [0.5, 0.5, 0.5], "si" => [0.25, 0.25, 0.25]]
    # written p before s: wannier90 orders by site, then l, then mr
    projections = parse_projections(["SI : p ; s"], lattice, atoms)
    @test length(projections) == 8
    @test [p.center for p in projections] == [fill([0.0, 0.0, 0.0], 4); fill([0.25, 0.25, 0.25], 4)]
    @test [(p.l, p.m) for p in projections[1:4]] == [(0, 1), (1, 1), (1, 2), (1, 3)]
    @test all(p -> p.n == 1 && p.α == 1 && p.zaxis == [0, 0, 1] && p.xaxis == [1, 0, 0], projections)

    # hybrids (l < 0) before s; a state listed twice counts once
    projections = parse_projections(["Ga: s; sp3; s"], lattice, atoms)
    @test [(p.l, p.m) for p in projections] == [(-3, 1), (-3, 2), (-3, 3), (-3, 4), (0, 1)]

    @test [(p.l, p.m) for p in parse_projections(["Ga:l=2,mr=5,1;l=-2"], lattice, atoms)] ==
        [(-2, 1), (-2, 2), (-2, 3), (2, 1), (2, 5)]
    @test [(p.l, p.m) for p in parse_projections(["Ga:dxy,fz(x2-y2),sp3d2-6"], lattice, atoms)] ==
        [(-5, 6), (2, 5), (3, 5)]
end

@testitem "parse_projections sites, units and modifiers" begin
    using LinearAlgebra
    lattice = [4.0 0 0; 0 5.0 0; 0 0 6.0]
    atoms = ["C" => [0.0, 0.0, 0.0]]
    p = only(parse_projections(["f=0.1,0.2,-0.3:s:r=2:zona=0.5"], lattice, atoms))
    @test p.center == [0.1, 0.2, -0.3]
    @test (p.n, p.α) == (2, 0.5)

    p = only(parse_projections(["c=1.0,2.5,3.0:s"], lattice, atoms))
    @test p.center ≈ [0.25, 0.5, 0.5]
    p = only(parse_projections(["Bohr", "c=1.0,2.5,3.0:s"], lattice, atoms))
    @test p.center ≈ [0.25, 0.5, 0.5] .* WannierIO.Bohr
    @test only(parse_projections(["ang", "c=1,0,0:s"], lattice, atoms)).center ≈ [0.25, 0, 0]

    # axes are normalized; a slightly nonorthogonal x is projected
    p = only(parse_projections(["C:px:z=0,0,2:x=1,0,0.001"], lattice, atoms))
    @test p.zaxis == [0, 0, 1]
    @test p.xaxis ≈ [1, 0, 0]
    @test abs(dot(p.zaxis, p.xaxis)) < 1.0e-12
    # an axially symmetric orbital gets an orthogonal x (wannier90: random)
    p = only(parse_projections(["C:pz:z=1,0,0"], lattice, atoms))
    @test p.zaxis == [1, 0, 0]
    @test abs(dot(p.zaxis, p.xaxis)) < 1.0e-12 && norm(p.xaxis) ≈ 1
    @test_throws ArgumentError parse_projections(["C:px:z=1,0,0"], lattice, atoms)

    @test isempty(parse_projections(String[], lattice, atoms))
end

@testitem "parse_projections spinors" begin
    lattice = [4.0 0 0; 0 4.0 0; 0 0 4.0]
    atoms = ["Fe" => [0.0, 0.0, 0.0]]
    projections = parse_projections(["Fe:s;pz"], lattice, atoms; spinors = true)
    @test projections isa Vector{WannierIO.SpinorHydrogenOrbital}
    @test [(p.l, p.spin) for p in projections] == [(0, 1), (0, -1), (1, 1), (1, -1)]
    @test all(p -> p.spin_qaxis == [0, 0, 1], projections)

    projections = parse_projections(["Fe:d:z=0,0,1(d)[1,0,0]"], lattice, atoms; spinors = true)
    @test length(projections) == 5
    @test all(p -> p.spin == -1 && p.spin_qaxis == [1, 0, 0], projections)
    @test [p.spin for p in parse_projections(["Fe:s(u)"], lattice, atoms; spinors = true)] == [1]
    # the parentheses of an f orbital name are not a spin selection
    @test length(parse_projections(["Fe:fz(x2-y2)"], lattice, atoms; spinors = true)) == 2
    @test length(parse_projections(["Fe:fz(x2-y2):z=0,0,1"], lattice, atoms)) == 1
end

@testitem "parse_projections errors" begin
    lattice = [4.0 0 0; 0 4.0 0; 0 0 4.0]
    atoms = ["Fe" => [0.0, 0.0, 0.0]]
    for lines in (
            ["random"], ["random", "Fe:s"], ["Fe"], ["Co:s"], ["Fe:q"], ["Fe:l=4"],
            ["Fe:l=1,mr=4"], ["Fe:l=1,m=1"], ["Fe:s(u)"], ["Fe:s[0,0,1]"],
            ["Fe:s:y=1,0,0"], ["f=0,0:s"],
        )
        @test_throws ArgumentError parse_projections(lines, lattice, atoms)
    end
    for lines in (["Fe:s(u):z=0,0,1"], ["Fe:s()"], ["Fe:s(u"], ["Fe:s[0,0,1](u)"])
        @test_throws ArgumentError parse_projections(lines, lattice, atoms; spinors = true)
    end
end

@testitem "parse_projections win" begin
    using LazyArtifacts
    win = read_win(artifact"Si2_valence/Si2_valence.win")
    nnkp = read_nnkp(artifact"Si2_valence/outputs/Si2_valence.nnkp")
    projections = parse_projections(win)
    @test length(projections) == length(nnkp.projections)
    for (p, q) in zip(projections, nnkp.projections)
        @test (p.n, p.l, p.m, p.α) == (q.n, q.l, q.m, q.α)
        @test p.center ≈ q.center atol = 1.0e-5
        @test p.zaxis ≈ q.zaxis && p.xaxis ≈ q.xaxis
    end

    @test_throws ArgumentError parse_projections(delete!(copy(win), "projections"))
    win["num_wann"] = length(projections) + 1
    @test_throws ArgumentError parse_projections(win)
end

@testitem "parse_projections spinor win" begin
    using LazyArtifacts
    win = read_win(artifact"Fe_soc/Fe.win")
    nnkp = read_nnkp(artifact"Fe_soc/outputs/Fe.nnkp")
    projections = parse_projections(win)
    @test projections isa Vector{WannierIO.SpinorHydrogenOrbital}
    @test length(projections) == length(nnkp.spinor_projections)
    for (p, q) in zip(projections, nnkp.spinor_projections)
        @test (p.n, p.l, p.m, p.α, p.spin) == (q.n, q.l, q.m, q.α, q.spin)
        @test p.center ≈ q.center atol = 1.0e-5
        @test p.spin_qaxis ≈ q.spin_qaxis
    end
end

@testitem "parse_projections crystal" begin
    using WannierIO: Crystal
    lattice = [0 1 1; 1 0 1; 1 1 0] * 2.715
    atoms = ["Si" => [0.0, 0.0, 0.0], "Ga" => [0.5, 0.5, 0.5], "Si" => [0.25, 0.25, 0.25]]
    crystal = Crystal(lattice, atoms)
    lines = ["SI:p;s", "Ga:sp3:z=1,1,0:x=1,-1,0"]
    @test parse_projections(lines, crystal) == parse_projections(lines, lattice, atoms)
    @test parse_projections(lines, crystal; spinors = true) ==
        parse_projections(lines, lattice, atoms; spinors = true)
    # `c=` goes through the lattice, a static matrix in the crystal
    lines = ["Bohr", "c=1.0,2.5,3.0:s"]
    @test only(parse_projections(lines, crystal)).center ≈
        only(parse_projections(lines, lattice, atoms)).center
end

@testitem "format_projections round trip" begin
    lattice = [4.0 0 0; 0 4.0 0; 0 0 4.0]
    inner(block) = split(block, "\n"; keepempty = false)[2:(end - 1)]
    projections = [
        WannierIO.HydrogenOrbital([0.1, 0.2, 0.3], 1, 2, 4, 1.0, [0, 0, 1], [1, 0, 0]),
        WannierIO.HydrogenOrbital([0.5, 0, 0], 2, -3, 2, 0.5, [0, 1, 0], [0, 0, 1]),
    ]
    block = format_projections(projections)
    @test startswith(block, "begin projections\n") && endswith(block, "end projections\n")
    @test occursin(":r=2:zona=0.5", block)
    @test parse_projections(inner(block), lattice, []) == projections

    spinor = [
        WannierIO.SpinorHydrogenOrbital([0, 0, 0], 1, 1, 1, 1.0, [0, 0, 1], [1, 0, 0], 1, [0, 0, 1]),
        WannierIO.SpinorHydrogenOrbital([0, 0, 0], 1, 1, 1, 1.0, [0, 0, 1], [1, 0, 0], -1, [0, 0, 1]),
        WannierIO.SpinorHydrogenOrbital([0, 0, 0], 1, 0, 1, 1.0, [0, 0, 1], [1, 0, 0], -1, [1, 0, 0]),
    ]
    block = format_projections(spinor)
    @test length(inner(block)) == 2
    @test occursin("(u,d)[", block) && occursin("(d)[", block)
    @test parse_projections(inner(block), lattice, []; spinors = true) == spinor
end
