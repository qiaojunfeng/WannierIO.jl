function istoml(io::IO)
    content = read(io, String)
    try
        TOML.parse(content)
    catch err
        err isa TOML.ParserError || rethrow()
        return false
    else
        return true
    end
end

function istoml(filename::AbstractString)
    return open(filename) do io
        istoml(io)
    end
end

"""
    $(SIGNATURES)

Convert types for TOML writing.

# Examples
```julia
using TOML
x = Dict(:a => (; b=1))
TOML.print(to_toml, x)
```
"""
function to_toml(x)
    x isa Pair && return Dict(x)
    # x isa Symbol && return String(x)
    x isa Mat3 && return eachcol(x)
    # use pairs other than Dict to preserve the order as much as possible
    x isa NamedTuple && return pairs(x)
    x isa HydrogenOrbital && return pairs(NamedTuple(x))
    # I explicitly throw an error here, to make sure that I don't forget to
    # handle a type. Sometimes, TOML.print just print a `false` for some type.
    return error("unhandled type $(typeof(x))")
end

"""
    $(SIGNATURES)

Write `kwargs` into `io` as a TOML file.

This is more convenient than `TOML.print`, in that:
1. Do some type conversion before writing.
2. Keep the order of kwargs, whereas `TOML.print` always requires a `Dict` which
    will change the order of keys.

# Examples
```julia
write_toml(stdout; b=2, a=1)
```
"""
function write_toml end

function write_toml(io, params::AbstractDict)
    return TOML.print(to_toml, io, params)
end

"""
I store atoms_frac and kpoint_path as Vector of StringVec3.
However, TOML.print does not accept Pair (specifically, StringVec3);
instead, I convert StringVec3 to Dict in write_win with toml format.
On reading I convert it back.
"""
function from_toml(d::AbstractDict)
    # StringVec3 are converted to Dict of length 1 when writing
    if length(d) == 1
        k, v = only(d)
        if isa(k, AbstractString) && isa(v, AbstractVector{<:Real})
            return stringvec3(k, v)
        end
    end

    # Need to do the conversion recursively
    for (k, v) in pairs(d)
        d[k] = from_toml(v)
    end
    return d
end

function from_toml(v::AbstractVector{<:Real})
    if length(v) == 3
        return vec3(v)
    else
        return v
    end
end

function from_toml(v::AbstractVector)
    # For a 3-vector of 3-vector, TOML returned type is Vector{Any}.
    # We want to convert it to Mat3.
    if length(v) == 3
        if all(x -> isa(x, AbstractVector{<:Real}) && length(x) == 3, v)
            return mat3(v...)
        end
    end
    return map(from_toml, v)
end

"""Fallback to doing nothing."""
from_toml(x) = x

function read_toml(io::IO)
    d = TOML.parse(read(io, String))
    d = from_toml(d)

    # TOML is unordered, we just return Dict.
    # Need to set value to Any, otherwise value type can be too narrow and then
    # I cannot assign Mat3 to it.
    return d::Dict{String,Any}
end
