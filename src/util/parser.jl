
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
        @assert !eof(io) "unexpected end of file"
        line = strip(readline(io))
        splitted = split(line)
        n_splitted = length(splitted)
        vec[(counter + 1):(counter + n_splitted)] = parse.(T, splitted)
        counter += n_splitted
    end

    return vec
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
```
"""
function parse_indices(str::AbstractString)
    segments = split(replace(strip(str), r"\s+" => ""), ",")
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
