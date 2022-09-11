"""
    get_recip_lattice(lattice::Mat3)

Return reciprocal lattice.
"""
get_recip_lattice(lattice::Mat3) = 2π * inv(lattice)'

"""
    get_lattice(recip_lattice::Mat3)

Return lattice.
"""
get_lattice(recip_lattice::Mat3) = inv(recip_lattice / (2π))'
