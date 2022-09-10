# Documentation

## Writing

The documentation is written in Markdown, then processed by `Documenter.jl`.

## Build

Build docs locally

```shell
julia --project=. make.jl; python -m http.server --directory build
```
