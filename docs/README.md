# Documentation

## Writing

The documentation is written in Markdown and built with
`Documenter.jl` + `DocumenterVitepress.jl`.

## Build

See [`make_serve.sh`](./make_serve.sh).

The static site is generated under `docs/build/1` for local builds.

Useful commands:

```bash
# from repository root
julia --project=docs docs/make.jl

# preview locally
cd docs
python -m http.server --directory build/1
```
