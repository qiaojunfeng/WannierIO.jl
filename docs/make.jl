using Documenter
using WannierIO

# Generate the HTML pages by Documenter.jl
makedocs(;
    sitename="WannierIO.jl",
    authors="Junfeng Qiao and contributors.",
    clean=true,
    modules=[WannierIO],
    pages=[
        "Home" => "index.md",
        "API" => [
            "Convention" => "api/convention.md",
            "Utilities" => "api/util.md",
            "Wannier90" => "api/w90.md",
            "Volumetric data" => "api/volumetric.md",
            "Quantum ESPRESSO" => "api/qe.md",
            "EPW" => "api/epw.md",
            "Index" => "api/index.md",
        ],
    ],
)

# Documenter will auto detect build environment; on local machine it will be
# skipped, so it's safe to run this script
deploydocs(; repo="github.com/qiaojunfeng/WannierIO.jl.git", devbranch="main")
