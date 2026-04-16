using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    # Compact in-memory fixtures that exercise common I/O and operator paths.
    win_txt = """
num_wann = 2
mp_grid = 1 1 1
begin unit_cell_cart
ang
5.43 0.0 0.0
0.0 5.43 0.0
0.0 0.0 5.43
end unit_cell_cart
begin atoms_frac
Si 0.0 0.0 0.0
Si 0.25 0.25 0.25
end atoms_frac
begin kpoints
0.0 0.0 0.0
end kpoints
begin projections
Si:s;p
end projections
begin kpoint_path
G  0.000 0.000 0.000    X  0.500 0.000 0.500
end kpoint_path
"""

    lattice = mat3(
        [5.43, 0.0, 0.0],
        [0.0, 5.43, 0.0],
        [0.0, 0.0, 5.43],
    )
    Rvectors = [vec3(0, 0, 0), vec3(1, 0, 0)]
    Rdegens = [1, 1]
    H = stack([
        ComplexF64[0.5 + 0.0im 0.1 + 0.0im; 0.1 + 0.0im 0.6 + 0.0im],
        ComplexF64[0.0 + 0.0im 0.2 + 0.0im; 0.2 + 0.0im 0.0 + 0.0im],
    ])
    rx = H
    ry = H
    rz = H

    hrdat = HrDat("precompile", Rvectors, Rdegens, H)
    tbdat = TbDat("precompile", lattice, Rvectors, Rdegens, H, rx, ry, rz)

    dense = stack([
        ComplexF64[1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 2.0 + 0.0im],
        ComplexF64[0.0 + 0.0im 3.0 + 0.0im; 3.0 + 0.0im 0.0 + 0.0im],
    ])

    cube = Cube(
        [vec3(0.0, 0.0, 0.0)],
        [14],
        vec3(0.0, 0.0, 0.0),
        mat3([0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5]),
        [0.0, 1.0],
        [0.0, 1.0],
        [0.0, 1.0],
        reshape(collect(1.0:8.0), 2, 2, 2),
    )

    xsf = Xsf(
        lattice,
        nothing,
        ["14"],
        [vec3(0.0, 0.0, 0.0)],
        vec3(0.0, 0.0, 0.0),
        lattice,
        [0.0, 1.0],
        [0.0, 1.0],
        [0.0, 1.0],
        reshape(collect(1.0:8.0), 2, 2, 2),
    )

    bxsf = Bxsf(
        0.0,
        vec3(0.0, 0.0, 0.0),
        lattice,
        [0.0, 1.0],
        [0.0, 1.0],
        [0.0, 1.0],
        reshape(collect(1.0:8.0), 1, 2, 2, 2),
    )

    serialized_text(writer) = begin
        io = IOBuffer()
        writer(io)
        String(take!(io))
    end

    @compile_workload begin
        mktemp() do winfile, io
            write(io, win_txt)
            close(io)

            win = read_win(winfile)
            write_win(winfile, win)
        end

        hrdat_txt = serialized_text() do io
            write_w90_hr_dat(io, hrdat)
        end
        read_w90_hr_dat(IOBuffer(hrdat_txt))

        tbdat_txt = serialized_text() do io
            write_w90_tb_dat(io, tbdat)
        end
        read_w90_tb_dat(IOBuffer(tbdat_txt))

        pack_hr = pack(hrdat)
        HrDat(pack_hr)
        wsvec = WsvecDat("precompile", Rvectors, 2)

        mktempdir() do tmpdir
            hrfile = joinpath(tmpdir, "precompile_hr.dat")
            write_w90_hr(hrfile, hrdat, wsvec)
            read_w90_hr(hrfile)

            tbfile = joinpath(tmpdir, "precompile_tb.dat")
            write_w90_tb(tbfile, tbdat, wsvec)
            read_w90_tb(tbfile)
        end

        sp = sparsify(dense)
        densify(sp)

        cube_txt = serialized_text() do io
            write_cube(io, cube)
        end
        read_cube(IOBuffer(cube_txt))

        xsf_txt = serialized_text() do io
            write_xsf(io, xsf)
        end
        read_xsf(IOBuffer(xsf_txt))

        bxsf_txt = serialized_text() do io
            write_bxsf(io, bxsf)
        end
        read_bxsf(IOBuffer(bxsf_txt))
    end
end
