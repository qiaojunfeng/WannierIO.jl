"""
    $(SIGNATURES)

Compare two (same-type) structs.
"""
function _isapprox(a::T, b::T) where {T}
    for f in propertynames(a)
        va = getfield(a, f)
        vb = getfield(b, f)

        if va isa String
            va == vb || return false
        elseif va isa Vector
            if eltype(va) isa BitVector
                all(va .== vb) || return false
            else
                all(va .≈ vb) || return false
            end
        else
            va ≈ vb || return false
        end
    end
    return true
end
