# Conventions

!!! tip

    In most cases, the units of the function arguments and returns are in angstrom unit for lattice,
    and fractional w.r.t lattice for atomic positions, etc.

!!! tip

    The following abbreviations are used throughout the code and documentation:
    * `n_bands` for number of bands
    * `n_wann` for number of WFs
    * `n_kpts` for number of kpoints
    * `n_bvecs` for number of b-vectors
    * `n_atoms` for number of atoms
    * `A` for `amn` matrices
    * `M` for `mmn` matrices
    * `E` for `eig` matrices
    * `S` for `spn` matrices

!!! note

    In most cases, for arrays we adopt the convention that `n_bands` is the first index,
    `n_wann` is the second index, and `n_kpts` is the third index.
    For example, `A` for the `amn` matrices is a 3D array of size `(n_bands, n_wann, n_kpts)`.
