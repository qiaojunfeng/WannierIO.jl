using YAML

@testset "read win" begin
    yaml_path = joinpath(@__DIR__, "test_data/win.yaml")
    test_data = YAML.load_file(yaml_path)

    win = read_win(joinpath(FIXTURE_PATH, "si2.win"))

    types_float = [
        AbstractFloat,
        AbstractVector{<:AbstractFloat},
        AbstractVector{<:AbstractVector{<:AbstractFloat}},
        AbstractArray{<:AbstractFloat},
    ]

    # convert types so that YAML can read/write arrays instead of strings
    preprocess_params(params; read=true) = begin
        params_dict = Dict(pairs(params))
        for (k, v) in pairs(params_dict)
            # no need to convert Float to vector
            v isa AbstractFloat && continue
            # convert Vector, Array to nested vectors
            if any(v isa T for T in types_float)
                if read
                    params_dict[k] = vec2mat(v)
                else
                    params_dict[k] = mat2vec(v)
                end
            end
        end
        return params_dict
    end

    # how to generate reference data:
    GENERATE_YAML = false
    if GENERATE_YAML
        # 1) convert arrays to nested vectors
        yaml_dict = preprocess_params(win; read=false)
        # 2) write to YAML
        YAML.write_file(yaml_path, yaml_dict)
    end

    # 3) upon loading, convert nested vectors to arrays
    test_data = preprocess_params(test_data)

    function test_kpoint_path()
        test_kpath = test_data["kpoint_path"]
        for (i, path) in enumerate(test_kpath)
            win_path = win.kpoint_path[i]
            k1, k2 = path
            # Dict to pair
            k1 = first(k1)
            k2 = first(k2)

            (Symbol(k1.first) != win_path[1].first) && return false
            (k1.second ≉ win_path[1].second) && return false
            (Symbol(k2.first) != win_path[2].first) && return false
            (k2.second ≉ win_path[2].second) && return false
        end
        return true
    end
    @test test_kpoint_path()

    for (k, v) in pairs(win)
        # already tested
        k == :kpoint_path && continue
        if any(v isa T for T in types_float)
            @test v ≈ test_data[string(k)]
            continue
        end
        @test v == test_data[string(k)]
    end
end

@testset "read/write win" begin
    win = read_win(joinpath(FIXTURE_PATH, "si2.win"))

    tmpfile = tempname(; cleanup=true)
    write_win(tmpfile; win...)
    win2 = read_win(tmpfile)

    for (k, v) in pairs(win)
        @test v == win2[k]
    end
end
