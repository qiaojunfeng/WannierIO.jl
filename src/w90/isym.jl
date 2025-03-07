export read_w90_isym #, write_w90_isym


mutable struct Symops 
    s
    ft 
    t_rev
    invs
end 

mutable struct LittleGroup
    h::Vector{Matrix{ComplexF64}}
    isym::Vector{Int64}
end 

"""
    $(SIGNATURES)

Read `prefix.isym`.

# Return
"""
function read_w90_isym(filename::AbstractString)
    # I use stream io to write mmn, so I should use plain julia `open`
    #println("I am starting")


    res = open("$filename") do io
        header_len = 60
        header = read(io, FString{header_len})

        # from FString to String
        header = strip(String(header))
        #println(header)
        # gfortran default integer size = 4
        # https://gcc.gnu.org/onlinedocs/gfortran/KIND-Type-Parameters.html
        line = split(readline(io))
        Tint = Int32
        n_symm, n_spin = parse.(Tint, line[1:2])
        #println(n_symm, " " , n_spin)
        
        
        #
        #read all symmetry operations
        #
        s        = [zeros(Float64, 3, 3) for _ in 1:n_symm]
        ft       = [zeros(Float64, 3) for _ in 1:n_symm]
        t_rev    = zeros(Int64, n_symm)
        invs     = zeros(Int64, n_symm) 

        for isym in 1:n_symm
            header = read(io, FString{header_len})
            #read sym[isym]
            for j in 1:3
                line = split(readline(io))
                s[isym][j,:] = parse.(Float64, line)
            end
            #println(s[isym])
            #read t[isym]
            line = split(readline(io))
            ft[isym]= parse.(Float64, line)
            #println(ft[isym])
            #read invs
            t_rev[isym]    = parse(Int64, readline(io))
            invs[isym] = parse(Int64, readline(io))
        end

        #
        #read k points in ibz
        #
        line = readline(io) #empty line
        line = readline(io) #line K points
        nkpoints = parse(Int64, readline(io))
        #println(nkpoints)
        kpoints = [zeros(Float64, 3) for _ in 1:nkpoints]
        for ik in 1:nkpoints
            kpoints[ik][:] = parse.(Float64, split(readline(io)))
            #println(kpoints[ik][:])
        end 

        #
        #read operations in the small group of each k point for symmetrization
        #
        line = readline(io) #empty line
        line = readline(io) #line Rep Mat .. 
        n_bands, n_blocks = parse.(Int64, split(readline(io)))
        #println(n_bands, n_blocks)

        gk = [LittleGroup(Vector{Vector{Matrix{ComplexF64}}}(), Vector{Vector{Int64}}()) for _ in 1:nkpoints]

        for ik in 1:nkpoints
            [gk[ik].h    = Vector{Matrix{ComplexF64}}() for _ in 1:nkpoints]
            [gk[ik].isym = Vector{Int64}() for _ in 1:nkpoints]
        end

        for iblocks in 1:n_blocks
            ik, isym, num_el  = parse.(Int64, split(readline(io)))
            push!(gk[ik].h,    zeros(ComplexF64, n_bands, n_bands))
            push!(gk[ik].isym, isym)
            for i_el in 1:num_el
                line = split(readline(io))
                m, n       = parse.(Int64, line[1:2])
                real, imag = parse.(Float64, line[3:4])
                gk[ik].h[end][m,n] = real + im * imag
            end
            #println(ik, " ", isym, " ", num_el)
        end
        #println("h read")
        #
        #read rotation matrix D_mn(k)
        #
        line = readline(io) #empty line
        line = readline(io) #line Rotation marix of Wannier functions
        line = split(readline(io))
        n_wann = parse(Int64, line[1])

        D = [zeros(ComplexF64, n_wann, n_wann) for _ in 1:n_symm]
        #println("Reading D")
        for isym in 1:n_symm
            isym_, n_blocks = parse.(Int64, split(readline(io)))
            @assert isym == isym_ "isym and isym in file do not correspond"
            for iblocks in 1:n_blocks
                line = split(readline(io))
                m, n       = parse.(Int64, line[1:2])
                real, imag = parse.(Float64, line[3:4])
                D[isym][m,n] = real + im * imag
            end
            #println(isym, " ", D[isym][:,:])
        end

        symops = []
        for isym in 1:n_symm
            push!(symops, Symops(s[isym], ft[isym], t_rev[isym], invs[isym]))
        end
    

        #println("Read D")
        @info "Reading isym file" filename n_symm nkpoints
        return (; symops, kpoints, gk, D)
    end


    return (res)
end
