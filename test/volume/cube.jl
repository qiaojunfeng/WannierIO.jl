@testitem "read/write cube" begin
    using LazyArtifacts
    const Bohr::Float64 = 0.52917721092 
    cube = read_cube(artifact"Si2_valence/reference/Si2_valence_00001.cube")
    @test cube.origin ≈ [-7.10461, -9.47281, -9.47281]*Bohr 
    @test cube.voxel_vectors ≈
        transpose([0.0 0.39470 0.39470; 0.39470 0.00000 0.39470; 0.39470 0.39470 0.00000])*Bohr 
    @test cube.X ≈ range(0, 19, 20)
    @test cube.Y ≈ cube.X
    @test cube.Z ≈ cube.X
    @test size(cube.W) == (20,20,20)
    @test size(cube.atom_positions) == (3,8)
    print(cube.atom_positions[:,1])
    @test cube.atom_positions[:,1] ≈ transpose([0.00000 -5.13111 -5.13111])*Bohr

    tmpfile = tempname(; cleanup=true)
    write_cube(tmpfile, cube.atom_positions, cube.atom_numbers, cube.origin, cube.voxel_vectors, cube.W)
    cube2 = read_cube(tmpfile)

    @test cube.atom_positions ≈ cube2.atom_positions
    @test cube.atom_numbers ≈ cube2.atom_numbers
    @test cube.origin ≈ cube2.origin
    @test cube.voxel_vectors ≈ cube2.voxel_vectors
    @test cube.X ≈ cube2.X
    @test cube.Y ≈ cube2.Y
    @test cube.Z ≈ cube2.Z
    @test cube.W ≈ cube2.W
end
