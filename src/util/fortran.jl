"""
    $(SIGNATURES)

Check if a sequence of chars is binary.
"""
function isbinary(chars::AbstractVector{UInt8})
    # normal ASCII chars
    text_chars = Vector{UInt8}([7, 8, 9, 10, 12, 13, 27])
    append!(text_chars, 0x20:0x99)
    # null character
    push!(text_chars, 0x00)
    # ASCII DEL (Delete) control character
    deleteat!(text_chars, text_chars .== 0x7F)

    # display([Char(_) for _ in chars])

    # Check without mutating the input.
    return any(x -> x ∉ text_chars, chars)
end

"""
    $(SIGNATURES)

Check if the IO is in binary format.

This function preserves the stream position for seekable IO.
"""
function isbinary(io::IO)
    seekable = applicable(seek, io, 0)
    seekable || throw(ArgumentError("`isbinary(io)` requires seekable IO"))

    pos = position(io)
    try
        raw_data = read(io, 1024)
        return isbinary(raw_data)
    finally
        seek(io, pos)
    end
end

"""
    $(SIGNATURES)

Check if the file is in binary format.
"""
function isbinary(filename::AbstractString)
    return open(filename) do io
        return isbinary(io)
    end
end
