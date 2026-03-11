"""
    $(SIGNATURES)

Strip comments (`!` or `#`) from a line.

# Arguments
- `line::AbstractString`: input line

# Keyword Arguments
- `spaces::Bool=false`: if true, preserve leading and trailing whitespace.
"""
function strip_comment(line::AbstractString; spaces::Bool=false)
    i = findfirst(r"!|#", line)
    res = isnothing(i) ? line : line[1:(i.start - 1)]
    return spaces ? res : strip(res)
end

@inline readstrip(io::IO) = strip(readline(io))

"""
    $(SIGNATURES)

Read a line from the input stream, ignoring empty lines (containing only whitespace).

Returns the first non-empty line found, or empty string if EOF is reached.

# Arguments
- `io::IO`: input stream

# Keyword Arguments
- `comment::Bool=false`: if true, keep comments (starting with `!` or `#`)
- `lower::Bool=false`: if true, convert the returned string to lowercase

# Returns
- A string containing the first non-empty line, with leading and trailing whitespace removed, or empty string if EOF is reached.
"""
function nextline(io::IO; comment::Bool=false, lower::Bool=false)
    while !eof(io)
        line = readstrip(io)
        if !comment
            line = strip_comment(line)
        end
        if !isempty(line)
            return lower ? lowercase(line) : line
        end
    end
    return ""
end

"""
    $(SIGNATURES)

Parse a string as `Float64`.

The is capable of parsing Fortran outputs, e.g. `1.0D-10`, to the ordinary `1e-10`.
"""
function parse_float(s::AbstractString)
    if occursin("*", s)
        return NaN
    else
        return parse(Float64, replace(lowercase(strip(s)), "d" => "e"))
    end
end

"""
    $(SIGNATURES)

Parse a string as `bool`.

This is capable of parsing Fortran outputs, e.g., `.true.`, `.false.`, `true`, `T`.
"""
function parse_bool(s::AbstractString)
    s = replace(lowercase(strip(s)), "." => "")[1]  # only 1st char
    return s == 't' || s == '1'
end

"""
    $(SIGNATURES)

Parse an integer as `bool`.

- `0`: `false`
- `1` or `-1`: `true`
"""
function parse_bool(i::Integer)
    return i != 0
end

"""
    $(SIGNATURES)

Parse a vector of `n_elements` elements of type `T` from `io`.

# Arguments
- `io`: input stream
- `T`: type of elements
- `n_elements::Int`: total number of elements

# Examples

Suppose a file `demo.txt` has the following content:
```
1  2  3  4  5  6  7  8  9  10
11 12 13 14 15 16 17 18 19 20
21 22 23
```

Then the following code parses the file and return a vector filled with 1 to 23:
```julia-repl
julia> vector = open("demo.txt") do io
    parse_vector(io, Int, 23)
end
```

The number of elements in each line can be different.
"""
function parse_vector(io::IO, T::Type, n_elements::Integer)
    vec = zeros(T, n_elements)

    counter = 0
    while counter < n_elements
        !eof(io) || error("unexpected end of file")
        line = strip(readline(io))
        splitted = split(line)
        n_splitted = length(splitted)
        vec[(counter + 1):(counter + n_splitted)] = parse.(T, splitted)
        counter += n_splitted
    end

    return vec
end

@inline function parse_vector(parts::AbstractVector{<:AbstractString}, T::Type=Float64)
    map(x -> parse(T, x), parts)
end
@inline function parse_vector(line::AbstractString, T::Type=Float64)
    parse_vector(split(line), T)
end

"""
    $(SIGNATURES)

Parse a string of comma-separated indices or range into a vector of integers.

E.g., the `exclude_bands` tag of `win` file.

# Examples

```julia-repl
julia> parse_indices("1-2, 5,8 -10")
6-element Vector{Int64}:
  1
  2
  5
  8
  9
 10
julia> parse_indices("1-2 5,8 -10")
```
"""
function parse_indices(str::AbstractString)
    s = strip(str)
    # Remove leading and trailing whitespaces beside `-` sign
    s = replace(s, r"\s+-" => "-")
    s = replace(s, r"-\s+" => "-")
    # Replace `,` with whitespace
    s = replace(s, "," => " ")
    # Remove consecutive whitespaces
    s = replace(s, r"\s+" => " ")
    # Split by whitespace
    segments = split(s)
    indices = Vector{Int}()
    for s in segments
        if isempty(s)
            continue
        elseif occursin(r"^\d+-\d+$", s)
            start, stop = parse.(Int, split(s, "-"))
            push!(indices, start:stop...)
        elseif occursin(r"^\d+$", s)
            push!(indices, parse(Int, s))
        else
            throw(ArgumentError("invalid index: $s while parsing $str"))
        end
    end
    return indices
end

"""
    $(SIGNATURES)

Convert a vector of integers to a string of comma-separated indices or range.

E.g., the `exclude_bands` tag of `win` file.

# Examples

```julia-repl
julia> format_indices([1, 2, 5, 8, 9, 10])
"1-2 5 8-10"
```
"""
function format_indices(indices::AbstractVector{T}) where {T<:Integer}
    groups = Vector{Vector{T}}()
    for (i, n) in enumerate(indices)
        if (i > 1) && (n == indices[i - 1] + 1)
            push!(groups[end], n)
        else
            push!(groups, [n])
        end
    end
    result = map(groups) do grp
        if length(grp) == 1
            return string(grp[1])
        else
            return string(grp[1], "-", grp[end])
        end
    end
    result = join(result, ", ")
    return result
end
