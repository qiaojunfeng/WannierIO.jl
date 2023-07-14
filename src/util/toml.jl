"""
Write `kwargs` into `io` as a TOML file.

Do some type conversion before writing.
"""
function write_toml(io; kwargs...)
    TOML.print(io, kwargs) do x
        x isa Pair && return Dict(x)
        x isa Symbol && return String(x)
        x isa Mat3 && return eachcol(x)
        error("unhandled type $(typeof(x))")
    end
end
