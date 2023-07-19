@testitem "read tb ws" begin
    using LazyArtifacts
    tbdat = read_w90_tbdat(artifact"Si2_valence/reference/ws/Si2_valence_tb.dat")

    @test length(tbdat.Rvectors) == 279
    @test tbdat.Rvectors[1] == [-4, 0, 2]
    @test length(tbdat.H) == 279
    H1 = ComplexF64[
        0.00080451304+1.6092791e-9im -0.00043090632-8.9535509e-9im -0.00013727432+1.7830477e-10im -0.00043090075+2.7695627e-9im
        -0.0004309027-2.2505901e-9im 0.00080451732+3.968764e-9im -0.00013726424+5.6728313e-9im -0.00043091084-2.4100158e-9im
        -0.00029850062+3.3998938e-9im -0.00029849652+1.2885406e-9im 0.00053380915+6.5106419e-10im -0.00029849523-4.1141943e-9im
        -0.00043089878+6.7332791e-9im -0.00043091377+5.6168704e-9im -0.0001372618-4.1114441e-9im 0.00080451507-6.2301744e-9im
    ]
    @test tbdat.H[1] ≈ H1
    @test length(tbdat.Rdegens) == 279
    @test tbdat.Rdegens[1] == 3
    @test tbdat.lattice ==
        [0.0 2.715265 2.715265; 2.715265 0.0 2.715265; 2.715265 2.715265 0.0]
    @test length(tbdat.r_x) == length(tbdat.r_y) == length(tbdat.r_z) == 279
    r_x1 = ComplexF64[
        3.6338163e-8+5.0525602e-8im -4.0150663e-5+1.506446e-8im 9.0607115e-5-5.4813038e-8im -0.0017405343-1.1882662e-8im
        0.0011077989+2.2693247e-8im -6.4595459e-9+5.4150924e-9im 0.00025422758-1.6273951e-8im 0.0011078169+2.8407141e-8im
        0.00012821382-2.969072e-8im 0.00013798637-2.1224108e-9im 5.4586044e-8+2.3435702e-8im 0.00012827449+8.9856419e-9im
        -0.0017405579-4.523964e-8im -4.007513e-5+3.8697593e-8im 9.0627952e-5+5.942338e-8im -1.6888038e-8+2.5373747e-9im
    ]
    r_y2 = ComplexF64[
        -2.5598094e-10-3.8026624e-9im 0.00074239939-1.8853816e-8im -0.00017865222+3.344151e-8im 0.00017582479+1.5301678e-8im
        0.00023456283+1.0709899e-8im -3.2666715e-8+4.7958915e-9im 0.00021987434-1.973943e-9im -1.3153397e-5+9.7076299e-9im
        1.319145e-5-1.4937437e-8im -0.00021977894+3.1806597e-9im -3.5269004e-8-1.3967952e-8im -0.00023454332+3.0129073e-8im
        -0.00017582461+1.1143559e-8im 0.00017861969+4.0119904e-8im -0.00074232978-1.8613896e-8im -2.2552721e-8-4.7774167e-9im
    ]
    r_z3 = ComplexF64[
        -2.7538102e-8-2.5754529e-9im 0.00017582603-4.8125766e-8im 0.00074233622-6.6380593e-9im -0.00017862626+1.7228737e-8im
        -0.00017585023-1.101898e-8im 6.7641871e-10-2.4185959e-8im 0.00017860646+3.7820915e-8im -0.00074235898-3.711759e-8im
        0.00023460039-2.0712049e-8im -1.318695e-5-4.598068e-9im 9.1295798e-9+4.1513854e-11im 0.00021978738+2.8018558e-8im
        1.3196588e-5+3.6553635e-9im -0.00023460902+3.2222708e-8im -0.00021979792-1.031378e-8im -3.1967074e-8-2.3895402e-8im
    ]
    @test tbdat.r_x[1] ≈ r_x1
    @test tbdat.r_y[2] ≈ r_y2
    @test tbdat.r_z[3] ≈ r_z3
end

@testitem "read tb mdrs" begin
    using LazyArtifacts
    # The two files are identical
    tbdat = read_w90_tbdat(artifact"Si2_valence/reference/mdrs/Si2_valence_tb.dat")
    tbdat_ws = read_w90_tbdat(artifact"Si2_valence/reference/ws/Si2_valence_tb.dat")

    for p in propertynames(tbdat)
        p == :header && continue
        @test tbdat_ws[p] ≈ tbdat[p]
    end
end
