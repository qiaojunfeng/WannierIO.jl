#!/bin/bash
# Build docs and start a local server to view them.
#
# There are two ways: LiveServer.jl or python's http.server
#   ./make_serve.sh
#   ./make_serve.sh py

USE_PYTHON=false
if [[ "$1" == "py" ]]; then
    USE_PYTHON=true
fi

# https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

if [[ $USE_PYTHON == false ]]; then
    # 1. Build docs then serve the Vitepress output directory.
    cd "$SCRIPT_DIR"
    julia --project make.jl
    # Use `0.0.0.0` to listen on all interfaces, so that port forward works.
    julia --project -e 'using LiveServer; serve(dir="build/1", host="0.0.0.0")'
else
    # 2. Use python's http.server to set up a local server.
    cd "$SCRIPT_DIR"
    julia --project make.jl
    python -m http.server --directory build/1
fi
