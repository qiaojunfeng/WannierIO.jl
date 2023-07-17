"""
    $(SIGNATURES)

Compute reciprocal lattice from lattice.
"""
get_recip_lattice(lattice::Mat3) = 2π * inv(lattice)'

"""
    $(SIGNATURES)

Compute lattice from reciprocal lattice.
"""
get_lattice(recip_lattice::Mat3) = inv(recip_lattice / (2π))'
