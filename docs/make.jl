using Documenter
using WannierIO

# Generate the HTML pages by Documenter.jl
makedocs(;
    sitename="WannierIO.jl",
    authors="Junfeng Qiao and contributors.",
    clean=true,
    modules=[WannierIO],
    pages=[
        "WannierIO" => [
            "Home" => "WannierIO/index.md",
            "API" => [
                "Convention" => "WannierIO/api/convention.md",
                "Util" => "WannierIO/api/util.md",
                "Wannier90" => "WannierIO/api/w90.md",
                "Volumetric data" => "WannierIO/api/volumetric.md",
                "QE" => "WannierIO/api/qe.md",
            ],
        ],
    ],
)
