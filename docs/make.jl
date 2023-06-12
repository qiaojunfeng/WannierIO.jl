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
        ],
    ],
)
