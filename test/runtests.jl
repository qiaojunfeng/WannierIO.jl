# using Aqua
# only test ambiguities in current module, otherwise fails due to ambiguities in other packages
#   https://github.com/JuliaTesting/Aqua.jl/issues/77
# disable project_extras since we don't use julia < 1.2
# Aqua.test_all(WannierIO; ambiguities=false, project_extras=false)
# Aqua.test_ambiguities(WannierIO)

using TestItemRunner

@run_package_tests verbose = true
