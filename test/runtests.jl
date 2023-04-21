using WannierIO
using Test

using Aqua
# only test ambiguities in current module, otherwise fails due to ambiguities in other packages
#   https://github.com/JuliaTesting/Aqua.jl/issues/77
# disable project_extras since we don't use julia < 1.2
Aqua.test_all(WannierIO; ambiguities=false, project_extras=false)
Aqua.test_ambiguities(WannierIO)

const TEST_PATH = @__DIR__
const FIXTURE_PATH = joinpath(TEST_PATH, "fixtures")

mat2vec(A::AbstractMatrix) = [v for v in eachcol(A)]

@testset "WannierIO.jl" begin
    include("util/fortran.jl")

    include("w90/win.jl")
    include("w90/wout.jl")
    include("w90/nnkp.jl")
    include("w90/amn.jl")
    include("w90/mmn.jl")
    include("w90/eig.jl")
    # include("w90/spn.jl")
    include("w90/unk.jl")
    include("w90/chk.jl")
    include("w90/band.jl")
    include("w90/tb.jl")

    include("qe/xml.jl")

    include("volume/bxsf.jl")
end
