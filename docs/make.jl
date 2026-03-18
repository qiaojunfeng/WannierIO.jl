using Documenter
using DocumenterVitepress
using WannierIO

# Generate docs with a Vitepress backend.
makedocs(;
    sitename="WannierIO.jl",
    authors="Junfeng Qiao and contributors.",
    clean=true,
    modules=[WannierIO],
    format=DocumenterVitepress.MarkdownVitepress(;
        repo="github.com/qiaojunfeng/WannierIO.jl", devbranch="main", devurl="dev"
    ),
    pages=[
        "Home" => "index.md",
        "Introduction" => [
            "Read/write files" => "intro/read_write.md",
            "Tight-binding operators" => "intro/tb.md",
            "Sparse storage" => "intro/sparse.md",
        ],
        "API" => [
            "Convention" => "api/convention.md",
            "Utilities" => "api/util.md",
            "Wannier90" => "api/w90.md",
            "Tight-binding" => "api/tb.md",
            "Volumetric data" => "api/volumetric.md",
            "EPW" => "api/epw.md",
        ],
    ],
)

# DocumenterVitepress handles deployments separately from Documenter.jl.
DocumenterVitepress.deploydocs(;
    repo="github.com/qiaojunfeng/WannierIO.jl.git",
    target=joinpath(@__DIR__, "build"),
    branch="gh-pages",
    devbranch="main",
    push_preview=true,
)
