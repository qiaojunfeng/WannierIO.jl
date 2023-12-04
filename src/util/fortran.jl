
"""
    $(SIGNATURES)

Check if a sequence of chars is binary.
"""
function isbinary(chars::AbstractVector{UInt8})::Bool
    # normal ASCII chars
    text_chars = Vector{UInt8}([7, 8, 9, 10, 12, 13, 27])
    append!(text_chars, 0x20:0x99)
    # null character
    push!(text_chars, 0x00)
    deleteat!(text_chars, text_chars .== 0x7F)

    # remove normal ASCII
    filter!(x -> x ∉ text_chars, chars)

    # display([Char(_) for _ in chars])
    return length(chars) > 0
end

"""
    $(SIGNATURES)

Check if the file is in binary format.
"""
function isbinary(filename::AbstractString)
    raw_data = zeros(UInt8, 1024)

    io = open(filename)
    readbytes!(io, raw_data)
    close(io)

    return isbinary(raw_data)
end
