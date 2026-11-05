@testitem "read/write orb" begin
    using LinearAlgebra
    dir = joinpath(@__DIR__, "orb_testfiles")

    # Si, external s, p projectors (d excluded), no spin
    orb = read_orb(joinpath(dir, "Si.orb"))
    @test !orb.spinors
    @test length(orb.orbitals) == 8
    o = orb.orbitals[2]
    @test o.atom == 1 && o.species == "Si" && o.radial == 2 && o.l == 1 && o.label == "m= 1"
    @test o.center ≈ [0.25, 1.25, 0.25]
    # pz of the ylmr2 harmonics with the i^l factor, in the w90 basis: -pz
    @test o.coefficients == reshape(ComplexF64[-1, 0, 0], 3, 1)
    @test orb.orbitals[5].center ≈ [0, 0, 0]
    @test any(startswith("Atomic projectors"), orb.header)

    # GaAs, spin-orbit |l j m_j> projectors of the UPF
    orb = read_orb(joinpath(dir, "GaAs_soc.orb"))
    @test orb.spinors
    @test length(orb.orbitals) == 16
    @test [o.species for o in orb.orbitals] == [fill("Ga", 8); fill("As", 8)]
    o = orb.orbitals[3]
    @test o.l == 1 && o.label == "j=0.5 m_j=-0.5" && size(o.coefficients) == (3, 2)
    @test norm(o.coefficients) ≈ 1
    # every projector is normalized and those of one site are orthonormal
    C = [vec(o.coefficients) for o in orb.orbitals[1:8] if o.l == 1]
    @test [dot(a, b) for a in C, b in C] ≈ I

    # noncollinear without spin-orbit: real harmonics times spin up or down
    orb = read_orb(joinpath(dir, "Si_noncollinear.orb"))
    @test orb.spinors
    @test orb.orbitals[2].label == "m= 1 s_z=-0.5"
    @test orb.orbitals[2].coefficients == reshape(ComplexF64[0, 1], 1, 2)
    @test orb.orbitals[1].radial == orb.orbitals[2].radial

    for name in ("Si.orb", "GaAs_soc.orb", "Si_noncollinear.orb")
        original = read_orb(joinpath(dir, name))
        tmpfile = tempname(; cleanup = true)
        write_orb(tmpfile, original)
        @test read_orb(tmpfile) == original
    end
end
