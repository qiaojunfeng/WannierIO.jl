using TestItemRunner
using WannierIO
using Documenter

# Use local artifacts, if they are changed locally.
# macro artifact_str(path)
#     joinpath(Sys.homedir(), "git/WannierDatasets/datasets", path)
# end

# Get filter string from command line arguments if provided
# Usage: julia --project test/runtests.jl "Name of Test"
filter_name = length(ARGS) > 0 ? ARGS[1] : nothing

if isnothing(filter_name)
    println("Running all tests...")

    @run_package_tests verbose = true

    DocMeta.setdocmeta!(WannierIO, :DocTestSetup, :(using WannierIO); recursive=true)
    doctest(
        WannierIO;
        # fix=true,  # update all the output in `jldoctest`
    )
else
    println("Running specific test: $filter_name")

    @run_package_tests verbose = true filter = ti -> endswith(ti.filename, filter_name)
end
