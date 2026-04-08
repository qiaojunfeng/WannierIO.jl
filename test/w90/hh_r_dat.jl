#= TODO test with some real files
@testitem "write HH_R" begin
    n_wann = 2
    n_rvecs = 3
    H = zeros(ComplexF64, n_wann, n_wann, n_rvecs)
    for ir in 1:n_rvecs
        for j in 1:n_wann
            for i in 1:n_wann
                H[i, j, ir] = complex(100 * ir + 10 * j + i, -(100 * ir + 10 * j + i))
            end
        end
    end
    R = Int[
        0 1 -1
        0 0 1
        0 0 0
    ]

    hhr = WannierIO.HHRDat(H, R, nothing, "HH_R test header")
    tmpfile = tempname(; cleanup=true)
    WannierIO.write_HH_R_dat(tmpfile, hhr)

    lines = readlines(tmpfile)
    @test lines[1] == "HH_R test header"
    @test lines[2] == "2"
    @test lines[3] == "3"
    @test length(lines) == 3 + n_wann * n_wann * n_rvecs

    first_row = split(strip(lines[4]))
    @test parse.(Int, first_row[1:5]) == [0, 0, 0, 1, 1]
    @test parse(Float64, first_row[6]) ≈ real(H[1, 1, 1])
    @test parse(Float64, first_row[7]) ≈ imag(H[1, 1, 1])

    last_row = split(strip(lines[end]))
    @test parse.(Int, last_row[1:5]) == [-1, 1, 0, 2, 2]
    @test parse(Float64, last_row[6]) ≈ real(H[2, 2, end])
    @test parse(Float64, last_row[7]) ≈ imag(H[2, 2, end])
end

@testitem "write HH_R ndegen" begin
    n_rvecs = 17
    H = reshape(complex.(1:n_rvecs, 0), 1, 1, n_rvecs)
    R = vcat(collect(1:n_rvecs)', zeros(Int, 1, n_rvecs), zeros(Int, 1, n_rvecs))
    N = collect(1:n_rvecs)

    hhr = WannierIO.HHRDat(H, R, N, WannierIO.default_header())
    tmpfile = tempname(; cleanup=true)
    WannierIO.write_HH_R_dat(tmpfile, hhr)

    ndegen_file = tmpfile * ".ndegen"
    @test isfile(ndegen_file)

    ndegen_lines = readlines(ndegen_file)
    @test length(ndegen_lines) == 2
    @test length(split(ndegen_lines[1])) == 15
    @test length(split(ndegen_lines[2])) == 2

    N_out = parse.(Int, vcat(split.(ndegen_lines)...))
    @test N_out == N
end

@testitem "write HH_R invalid input" begin
    H = zeros(ComplexF64, 2, 3, 1)
    R = zeros(Int, 3, 1)
    hhr = WannierIO.HHRDat(H, R, nothing, "")
    @test_throws ArgumentError WannierIO.write_HH_R_dat(tempname(; cleanup=true), hhr)

    H2 = zeros(ComplexF64, 2, 2, 2)
    R2 = zeros(Int, 3, 2)
    N2 = [1]
    hhr2 = WannierIO.HHRDat(H2, R2, N2, "")
    @test_throws ArgumentError WannierIO.write_HH_R_dat(tempname(; cleanup=true), hhr2)
end
=#
