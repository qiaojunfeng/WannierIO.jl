using TestItemRunner

# Filter tests
# @run_package_tests verbose=true filter=ti->endswith(ti.filename, "type.jl")

@run_package_tests verbose=true
