using TestItemRunner

# Filter tests
# @run_package_tests verbose=true filter=ti->endswith(ti.filename, "type.jl")

# Use local artifacts, if they are changed locally.
# macro artifact_str(path)
#     joinpath(Sys.homedir(), "git/WannierDatasets/datasets", path)
# end

@run_package_tests verbose=true
