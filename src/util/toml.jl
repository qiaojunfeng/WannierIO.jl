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

function write_toml(io, params::NamedTuple)
    # Note that we cannot pass NamedTuple to TOML.print, but we can pass
    # pairs(::NamedTuple) which is a subtype of AbstractDict.
    return write_toml(io, pairs(params))
end

function write_toml(io; kwargs...)
    # typeof(kwargs) = pairs(::NamedTuple), and this keeps the order.
    return write_toml(io, kwargs)
end
