using Documenter
using WannierIO

makedocs(;
    sitename="WannierIO.jl",
    authors="Junfeng Qiao and contributors.",
    modules=[WannierIO],
    pages=[
        "Home" => "index.md",
        "API" => [
            "Convention" => "api/convention.md",
            "Util" => "api/util.md",
            "Wannier90" => "api/w90.md",
            "Volumetric data" => "api/volumetric.md",
            "QE" => "api/qe.md",
        ],
    ],
)
