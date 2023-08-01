@testitem "read r WS" begin
    using LazyArtifacts
    rdat = read_w90_rdat(artifact"Si2_valence/reference/WS/Si2_valence_r.dat")
    Rvectors1 = [-4, 0, 2]
    @test rdat.Rvectors[1] == Rvectors1
    r_x_end = ComplexF64[
        0.0-0.0im -0.001846-0.0im 0.000296+0.0im 0.002093+0.0im
        -0.000806-0.0im -0.0-0.0im -0.000149+0.0im -0.000806+0.0im
        -4.4e-5-0.0im 0.00034+0.0im 0.0-0.0im -4.4e-5+0.0im
        0.002093+0.0im -0.001846+0.0im 0.000296-0.0im -0.0-0.0im
    ]
    @test rdat.r_x[end] ≈ r_x_end

    r_y_end = ComplexF64[
        -0.0+0.0im -0.002093-0.0im -0.000296-0.0im 0.001846-0.0im
        -0.002093+0.0im -0.0-0.0im -0.000296-0.0im 0.001846+0.0im
        4.4e-5-0.0im 4.4e-5+0.0im -0.0+0.0im -0.00034+0.0im
        0.000806+0.0im 0.000806+0.0im 0.000149+0.0im -0.0-0.0im
    ]
    @test rdat.r_y[end] ≈ r_y_end

    r_z_end = ComplexF64[
        0.0-0.0im -0.000806-0.0im -0.000149+0.0im -0.000806+0.0im
        -0.001846-0.0im 0.0+0.0im 0.000296-0.0im 0.002093+0.0im
        0.00034+0.0im -4.4e-5+0.0im 0.0-0.0im -4.4e-5-0.0im
        -0.001846+0.0im 0.002093+0.0im 0.000296+0.0im -0.0+0.0im
    ]
    @test rdat.r_z[end] ≈ r_z_end
end

@testitem "read r MDRS" begin
    using LazyArtifacts
    # The two files are identical
    rdat = read_w90_rdat(artifact"Si2_valence/reference/MDRS/Si2_valence_r.dat")
    rdat_ws = read_w90_rdat(artifact"Si2_valence/reference/WS/Si2_valence_r.dat")

    for p in propertynames(rdat)
        p == :header && continue
        @test rdat_ws[p] ≈ rdat[p]
    end
end
