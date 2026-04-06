using TestItemRunner
using WannierIO
using Documenter

#= Use local artifacts, if they are changed locally.
macro artifact_str(path)
    joinpath(Sys.homedir(), "git/WannierDatasets/datasets", path)
end
=#

# Get filter strings from command line arguments if provided
# Usage: julia --project test/runtests.jl "w90/amn.jl" "util/parser.jl"
filter_names = isempty(ARGS) ? nothing : ARGS

if isnothing(filter_names)
    println("Running all tests...")

    @run_package_tests verbose = true

    DocMeta.setdocmeta!(WannierIO, :DocTestSetup, :(using WannierIO); recursive = true)
    doctest(
        WannierIO
        # fix=true,  # update all the output in `jldoctest`
    )
else
    println("Running specific tests: $(join(filter_names, ", "))")

    @run_package_tests verbose = true filter =
        ti -> any(name -> endswith(ti.filename, name), filter_names)
end
