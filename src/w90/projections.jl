export parse_projections, format_projections

# The `projections` block of a `win` file is kept as raw lines by `read_win`;
# `parse_projections` expands it into the orbitals that `wannier90.x -pp`
# writes to the `nnkp` file, following `w90_readwrite_get_projections` of
# wannier90 (src/readwrite.F90). `format_projections` is its inverse for
# explicit (`f=…:l=…,mr=…`) orbitals.

"""Number of `mr` states of the angular momentum `l` (`l < 0`: w90 hybrids)."""
_n_mr(l::Integer) = l >= 0 ? 2l + 1 : 1 - l

"""The `(l, mr)` states of each named orbital of the w90 projection syntax."""
const _W90_ORBITAL_STATES = let
    states = Dict{String, Vector{Tuple{Int, Int}}}(
        "s" => [(0, 1)],
        "p" => [(1, m) for m in 1:3],
        "pz" => [(1, 1)], "px" => [(1, 2)], "py" => [(1, 3)],
        "d" => [(2, m) for m in 1:5],
        "dz2" => [(2, 1)], "dxz" => [(2, 2)], "dyz" => [(2, 3)],
        "dx2-y2" => [(2, 4)], "dxy" => [(2, 5)],
        "f" => [(3, m) for m in 1:7],
        "fz3" => [(3, 1)], "fxz2" => [(3, 2)], "fyz2" => [(3, 3)], "fxyz" => [(3, 4)],
        "fz(x2-y2)" => [(3, 5)], "fx(x2-3y2)" => [(3, 6)], "fy(3x2-y2)" => [(3, 7)],
    )
    for (l, name) in zip(-1:-1:-5, ("sp", "sp2", "sp3", "sp3d", "sp3d2"))
        states[name] = [(l, m) for m in 1:_n_mr(l)]
        for m in 1:_n_mr(l)
            states["$name-$m"] = [(l, m)]
        end
    end
    states
end

"""Parse a comma-separated triplet `x,y,z`, as w90 `utility_string_to_coord`."""
function _parse_projection_coord(s::AbstractString)
    parts = split(s, ',')
    length(parts) == 3 ||
        throw(ArgumentError("expected three comma-separated numbers, got `$s`"))
    return Vec3(parse_float.(parts))
end

"""Fractional centers of the site part (before the first `:`) of a projection line."""
function _projection_sites(site::AbstractString, lattice, atoms, to_angstrom::Real)
    if startswith(site, "c=")
        cart = _parse_projection_coord(site[3:end]) .* to_angstrom
        return [Vec3(lattice \ cart)]
    elseif startswith(site, "f=")
        return [_parse_projection_coord(site[3:end])]
    end
    centers = [Vec3(frac) for (label, frac) in atoms if lowercase(label) == site]
    isempty(centers) && throw(
        ArgumentError("projection site `$site` is neither `c=`, `f=` nor an atom label")
    )
    return centers
end

"""
Remove a trailing `[…]` (open = '[') or `(…)` (open = '(') group from `s`.
Returns `(rest, content)`, `content === nothing` if there is no group.
"""
function _strip_trailing_group(s::AbstractString, open::Char, close::Char)
    i = findlast(open, s)
    isnothing(i) && return s, nothing
    j = findnext(close, s, i)
    isnothing(j) && throw(ArgumentError("no closing `$close` in projection `$s`"))
    j == lastindex(s) || throw(
        ArgumentError("unexpected text after `$close` in projection `$s`")
    )
    return s[1:prevind(s, i)], s[nextind(s, i):prevind(s, j)]
end

"""The set of `(l, mr)` states of the angular part of a projection line."""
function _projection_states(angular::AbstractString)
    states = Set{Tuple{Int, Int}}()
    for group in split(angular, ';')
        if startswith(group, "l=")
            parts = split(group, ',')
            l = parse(Int, parts[1][3:end])
            -5 <= l <= 3 || throw(ArgumentError("l = $l out of range -5:3 in `$group`"))
            if length(parts) == 1
                union!(states, (l, m) for m in 1:_n_mr(l))
                continue
            end
            startswith(parts[2], "mr=") ||
                throw(ArgumentError("expected `mr=` after `l=` in `$group`"))
            for s in [parts[2][4:end]; parts[3:end]]
                m = parse(Int, s)
                1 <= m <= _n_mr(l) ||
                    throw(ArgumentError("mr = $m out of range 1:$(_n_mr(l)) for l = $l"))
                push!(states, (l, m))
            end
        else
            for name in split(group, ',')
                haskey(_W90_ORBITAL_STATES, name) ||
                    throw(ArgumentError("unknown orbital `$name` in projection"))
                union!(states, _W90_ORBITAL_STATES[name])
            end
        end
    end
    return states
end

"""
Normalize the `zaxis` and `xaxis` of the orbital `(l, mr)` and make them
orthogonal, as w90 does: a nonorthogonal `xaxis` is projected onto the plane
normal to `zaxis` if they are within 1e-2 of orthogonal, and replaced for the
axially symmetric `pz`, `dz2`, `fz3` (where w90 draws a random one, here the
first of x̂, ŷ orthogonalized against `zaxis`).
"""
function _orthonormal_axes(zaxis, xaxis, l::Integer, m::Integer)
    z = zaxis / norm(zaxis)
    x = xaxis / norm(xaxis)
    c = dot(z, x)
    abs(c) > 1.0e-6 || return z, x
    if l >= 0 && m == 1
        for e in (Vec3(1.0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0))
            v = e - dot(e, z) * z
            norm(v) > 0.1 && return z, v / norm(v)
        end
    end
    abs(c) > 1.0e-2 && throw(
        ArgumentError("projection zaxis $zaxis and xaxis $xaxis are not orthogonal")
    )
    return z, (x - c * z) / sqrt(1 - c^2)
end

"""
    parse_projections(win)
    parse_projections(lines, lattice, atoms; spinors = false)

Expand the `projections` block of a `win` file into the orbitals that
`wannier90.x -pp` writes to the `nnkp` file: a `Vector{HydrogenOrbital}`, or
with `spinors` a `Vector{SpinorHydrogenOrbital}`.

The first form takes a parsed `win` ([`read_win`](@ref)) and checks that the
block defines at least `num_wann` orbitals. The second takes the raw lines of
the block, the lattice vectors as columns (Å), and the `atoms_frac` pairs
`label => fractional position`.

The syntax is that of wannier90, case insensitive and ignoring spaces: an
optional first line `Ang` or `Bohr` (the unit of `c=`), then one
`site:orbitals[:modifiers][(spins)][[spin axis]]` per line, where
- `site` is `c=x,y,z` (Cartesian), `f=x,y,z` (fractional), or an atom label,
    which expands to every atom of that label, in the order of `atoms`;
- `orbitals` are `;`-separated groups, each either `l=…` with an optional
    `,mr=…,…` list, or `,`-separated names (`s`, `p`, `px`, `dxy`, `fz3`,
    `sp3`, `sp3-1`, …);
- `modifiers` are `z=x,y,z`, `x=x,y,z`, `r=n` and `zona=α`, separated by `:`;
- `(u)`, `(d)` or `(u,d)` (default) selects the spin components and `[x,y,z]`
    the spin quantization axis (default `[0,0,1]`), both only with `spinors`.

As in wannier90, the orbitals of a line are ordered by site, then `l`
(hybrids `l < 0` first), then `mr`, then spin (up first), whatever the order
in which they are written; the axes are normalized and made orthogonal. A
nonorthogonal `x=` of the axially symmetric `pz`, `dz2`, `fz3` is replaced by
the first of x̂, ŷ orthogonalized against `z=`, where wannier90 draws a random
axis. Throws an `ArgumentError` for `random` projections, which wannier90
draws from an unseeded generator, and, more strictly than wannier90, for
text after a closing `)` or `]` and for unknown modifiers, which it ignores.

See also [`format_projections`](@ref) for the inverse.
"""
function parse_projections end

function parse_projections(win::AbstractDict)
    haskey(win, "projections") || throw(ArgumentError("the win has no projections block"))
    projections = parse_projections(
        win["projections"], win["unit_cell_cart"], win["atoms_frac"];
        spinors = get(win, "spinors", false),
    )
    length(projections) >= win["num_wann"] || throw(
        ArgumentError(
            "the projections block defines $(length(projections)) orbitals, " *
                "fewer than num_wann = $(win["num_wann"])"
        ),
    )
    return projections
end

function parse_projections(
        lines::AbstractVector{<:AbstractString}, lattice::AbstractMatrix,
        atoms::AbstractVector; spinors::Bool = false,
    )
    lines = filter(!isempty, [replace(lowercase(line), r"\s" => "") for line in lines])
    to_angstrom = 1.0
    if !isempty(lines) && !occursin(':', first(lines))
        unit = popfirst!(lines)
        if startswith(unit, "random")
            throw(ArgumentError("random projections are not supported"))
        elseif startswith(unit, "b")
            to_angstrom = Bohr
        elseif !startswith(unit, "a")
            throw(ArgumentError("malformed projection `$unit`"))
        end
    end
    T = spinors ? SpinorHydrogenOrbital : HydrogenOrbital
    projections = T[]
    for line in lines
        _parse_projection_line!(projections, line, lattice, atoms, to_angstrom)
    end
    return projections
end

"""Append the orbitals of one projection line to `projections`."""
function _parse_projection_line!(
        projections::AbstractVector{T}, line::AbstractString, lattice, atoms,
        to_angstrom::Real,
    ) where {T <: Orbital}
    spinors = T === SpinorHydrogenOrbital
    occursin(':', line) || throw(ArgumentError("malformed projection `$line`"))
    site, rest = split(line, ':'; limit = 2)
    centers = _projection_sites(site, lattice, atoms, to_angstrom)

    rest, qaxis = _strip_trailing_group(rest, '[', ']')
    isnothing(qaxis) || spinors ||
        throw(ArgumentError("spin quantization axis in `$line` but spinors = false"))
    spin_qaxis = isnothing(qaxis) ? Vec3(0.0, 0.0, 1.0) : _parse_projection_coord(qaxis)

    # the last `(…)` is a spin selection, unless it belongs to an f orbital name
    i = findlast('(', rest)
    is_forbital = !isnothing(i) &&
        any(startswith(rest[i:end], s) for s in ("(x2-y2)", "(x2-3y2)", "(3x2-y2)"))
    rest, spin = is_forbital ? (rest, nothing) : _strip_trailing_group(rest, '(', ')')
    isnothing(spin) || spinors ||
        throw(ArgumentError("spin selection in `$line` but spinors = false"))
    spins = if isnothing(spin)
        spinors ? (1, -1) : (0,)
    else
        up, dn = 'u' in spin, 'd' in spin
        up || dn || throw(ArgumentError("spin selection `($spin)` has neither u nor d"))
        up && dn ? (1, -1) : up ? (1,) : (-1,)
    end

    angular, modifiers = occursin(':', rest) ? split(rest, ':'; limit = 2) : (rest, "")
    states = sort!(collect(_projection_states(angular)))

    zaxis, xaxis = Vec3(0.0, 0.0, 1.0), Vec3(1.0, 0.0, 0.0)
    n, α = 1, 1.0
    for modifier in split(modifiers, ':'; keepempty = false)
        key, value = occursin('=', modifier) ? split(modifier, '='; limit = 2) : (modifier, "")
        if key == "z"
            zaxis = _parse_projection_coord(value)
        elseif key == "x"
            xaxis = _parse_projection_coord(value)
        elseif key == "r"
            n = parse(Int, value)
        elseif key == "zona"
            α = parse_float(value)
        else
            throw(ArgumentError("unknown projection modifier `$modifier` in `$line`"))
        end
    end

    for center in centers, (l, m) in states
        z, x = _orthonormal_axes(zaxis, xaxis, l, m)
        for spin in spins
            if spinors
                push!(projections, T(center, n, l, m, α, z, x, spin, spin_qaxis))
            else
                push!(projections, T(center, n, l, m, α, z, x))
            end
        end
    end
    return projections
end

_format_vector(v) = join([@sprintf("%.10f", abs(x) < 1.0e-12 ? 0.0 : x) for x in v], ",")

function _format_orbital(o::Orbital)
    s = "f=$(_format_vector(o.center)):l=$(o.l),mr=$(o.m)"
    s *= ":z=$(_format_vector(o.zaxis)):x=$(_format_vector(o.xaxis))"
    o.n == 1 || (s *= ":r=$(o.n)")
    o.α == 1 || (s *= ":zona=$(o.α)")
    return s
end

"""
    format_projections(projections) -> String

The `begin projections ... end projections` block of a `.win` file for the
`HydrogenOrbital`s `projections`, one `f=…:l=…,mr=…:z=…:x=…` line per
orbital (fractional centers, Cartesian axes; `:r=…` and `:zona=…` when not
the default 1). For `SpinorHydrogenOrbital`s (which need `spinors = .true.`),
consecutive spin-up and spin-down copies of one orbital with a common
quantization axis are written as one line with `(u,d)[…]`, single ones with
`(u)[…]` or `(d)[…]`.

See also [`parse_projections`](@ref) for the inverse.
"""
function format_projections(projections::AbstractVector{HydrogenOrbital})
    lines = [_format_orbital(o) for o in projections]
    return join(["begin projections"; lines; "end projections"], "\n") * "\n"
end

function format_projections(projections::AbstractVector{SpinorHydrogenOrbital})
    lines = String[]
    i = 1
    while i <= length(projections)
        o = projections[i]
        spins = o.spin > 0 ? "u" : "d"
        if i < length(projections)
            n = projections[i + 1]
            same = (n.center, n.n, n.l, n.m, n.α, n.zaxis, n.xaxis, n.spin_qaxis) ==
                (o.center, o.n, o.l, o.m, o.α, o.zaxis, o.xaxis, o.spin_qaxis)
            if same && o.spin > 0 && n.spin < 0
                spins = "u,d"
                i += 1
            end
        end
        push!(lines, "$(_format_orbital(o))($spins)[$(_format_vector(o.spin_qaxis))]")
        i += 1
    end
    return join(["begin projections"; lines; "end projections"], "\n") * "\n"
end
