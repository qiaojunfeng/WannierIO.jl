
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
