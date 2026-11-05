export AtomicOrbital, Orb, read_orb, write_orb

# The `prefix.orb` file is written by Quantum ESPRESSO's pw2wannier90.x when
# `atom_proj = .true.`: it lists the atomic projectors of the columns of the
# `.amn` (or `.iamn`) file. Each projector is
#   φ(r) = R(|r - τ|) Σ_{mr,s} c[mr, s] Y_{l,mr}(r - τ) |s⟩,
# with Y_{l,mr} the real spherical harmonics of the Wannier90 projection table
# in the global Cartesian axes and |s⟩ (s = 1, 2) spin up and down along z.
# R is real and shared by the projectors of one atom with the same radial key.
#
# Layout: comment lines starting with `#`, then `n_proj npol`, then per
# projector one line `ip atom species τ₁ τ₂ τ₃ radial_key l label…` (τ in
# fractional coordinates) followed by `2l + 1` lines of `npol` complex numbers
# `Re c[mr, s] Im c[mr, s]`.

"""
An atomic projector of pw2wannier90's `atom_proj`, as listed in a `prefix.orb`
file: a radial function times an exact combination of real spherical
harmonics and spin states,
```math
φ(\\mathbf{r}) = R(|\\mathbf{r} - \\boldsymbol{τ}|)
\\sum_{m_r, s} c_{m_r s} Y_{l m_r}(\\mathbf{r} - \\boldsymbol{τ}) |s⟩,
```
with ``Y_{l m_r}`` the real harmonics of the Wannier90 projection table in the
global Cartesian axes (`mr` as in [`HydrogenOrbital`](@ref)) and ``|s⟩`` spin
up (`s = 1`) and down (`s = 2`) along ``z``. The radial function is not
stored: projectors of the same `species` and `radial` key share it. For a
spin-orbit calculation the coefficients are those of the ``|l, j, m_j⟩``
spinors of the pseudopotential, which have no `HydrogenOrbital` counterpart.

$(TYPEDEF)

# Fields

$(FIELDS)
"""
@kwdef struct AtomicOrbital <: Orbital
    """3 real numbers of the projection center (the atom), in fractional coordinates"""
    center::Vec3{Float64}
    """index of the atom in the DFT calculation"""
    atom::Int
    """atomic species label"""
    species::String
    """radial function key: the same for the projectors of one species that
    share the radial function"""
    radial::Int
    """angular momentum ``l \\ge 0``"""
    l::Int
    """magnetic quantum numbers as printed by pw2wannier90 (`m= 1`,
    `m= 1 s_z= 0.5`, or `j=0.5 m_j=-0.5`)"""
    label::String
    """`(2l + 1) × n_spin` matrix ``c_{m_r s}``, `n_spin = 1` without and `2`
    with spinors"""
    coefficients::Matrix{ComplexF64}
end

function Base.:(==)(a::AtomicOrbital, b::AtomicOrbital)
    return a.center == b.center && a.atom == b.atom && a.species == b.species &&
        a.radial == b.radial && a.l == b.l && a.label == b.label &&
        a.coefficients == b.coefficients
end

"""
Container for a `prefix.orb` file of pw2wannier90.

$(TYPEDEF)

# Fields

$(FIELDS)
"""
struct Orb
    """the comment lines of the file, without the leading `#`"""
    header::Vector{String}
    """whether the projectors are spinors (two spin components)"""
    spinors::Bool
    """the projectors, in the order of the `amn` columns"""
    orbitals::Vector{AtomicOrbital}
end

Base.:(==)(a::Orb, b::Orb) = a.spinors == b.spinors && a.orbitals == b.orbitals

"""
    read_orb(file)

Read a `prefix.orb` file written by pw2wannier90.x with `atom_proj = .true.`:
the atomic projectors of the `amn` columns, see [`AtomicOrbital`](@ref).

# Arguments
- `file`: the file name, or an `IO`.

# Return
- an [`Orb`](@ref)
"""
function read_orb end

function read_orb(io::IO)
    header = String[]
    line = readline(io)
    while startswith(line, "#")
        push!(header, strip(line[2:end]))
        line = readline(io)
    end
    n_proj, npol = parse.(Int, split(line)[1:2])
    npol in (1, 2) || error("npol = $npol must be 1 or 2")
    orbitals = map(1:n_proj) do ip
        t = split(readline(io))
        length(t) >= 8 || error("projector $ip: expected `ip atom species τ₁ τ₂ τ₃ radial l label`")
        parse(Int, t[1]) == ip || error("projector $ip is numbered $(t[1])")
        l = parse(Int, t[8])
        l >= 0 || error("projector $ip: l = $l must be non-negative")
        c = zeros(ComplexF64, 2l + 1, npol)
        for mr in 1:(2l + 1)
            x = parse_float.(split(readline(io)))
            length(x) == 2npol || error("projector $ip, mr = $mr: expected $(2npol) numbers")
            c[mr, :] = x[1:2:end] .+ im .* x[2:2:end]
        end
        AtomicOrbital(;
            center = Vec3(parse_float.(t[4:6])), atom = parse(Int, t[2]), species = String(t[3]),
            radial = parse(Int, t[7]), l, label = join(t[9:end], " "), coefficients = c,
        )
    end
    return Orb(header, npol == 2, orbitals)
end

read_orb(filename::AbstractString) = open(read_orb, filename)

"""
    write_orb(file, orb)

Write an [`Orb`](@ref) in the `prefix.orb` format of pw2wannier90.

# Arguments
- `file`: the file name, or an `IO`.
- `orb`: the projectors.
"""
function write_orb end

function write_orb(io::IO, orb::Orb)
    for line in orb.header
        println(io, "# ", line)
    end
    npol = orb.spinors ? 2 : 1
    @printf(io, "%8d%8d\n", length(orb.orbitals), npol)
    for (ip, o) in enumerate(orb.orbitals)
        size(o.coefficients) == (2o.l + 1, npol) ||
            error("projector $ip: coefficients must be $(2o.l + 1) × $npol")
        @printf(io, "%6d%6d %3s", ip, o.atom, o.species)
        for x in o.center
            @printf(io, "%20.15f", x)
        end
        @printf(io, "%6d%6d  %s\n", o.radial, o.l, o.label)
        for mr in 1:(2o.l + 1)
            for s in 1:npol
                c = o.coefficients[mr, s]
                @printf(io, "%24.15e%24.15e", real(c), imag(c))
            end
            println(io)
        end
    end
    return nothing
end

write_orb(filename::AbstractString, orb::Orb) = open(io -> write_orb(io, orb), filename, "w")
