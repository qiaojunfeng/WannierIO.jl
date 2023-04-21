using EzXML

"""
    read_qe_xml(filename)

Read atomic structure and band structure from QE's XML output.

# Return
- `lattice`: `3 * 3`, Å, each column is a lattice vector
- `atom_positions`: `3 * n_atoms`, fractional, each column is a position
- `atom_labels`: `n_atoms`, each element is the label of the corresponding atom
- `recip_lattice`: `3 * 3`, Å⁻¹, each column is a reciprocal lattice vector
- `kpoints`: `3 * n_kpts`, fractional, each column is a kpoint
- `E`: `n_bands * n_kpts`, eV
- `fermi_energy`: eV
"""
function read_qe_xml(filename::AbstractString)
    # from qe/Modules/constants.f90
    BOHR_RADIUS_ANGS = 0.529177210903  # Angstrom
    HARTREE_SI = 4.3597447222071e-18 # J
    ELECTRONVOLT_SI = 1.602176634e-19     # J
    AUTOEV = HARTREE_SI / ELECTRONVOLT_SI

    doc = readxml(filename)
    output = findfirst("/qes:espresso/output", root(doc))

    # atoms
    atomic_structure = findfirst("atomic_structure", output)
    alat = parse(Float64, atomic_structure["alat"])
    # from bohr to angstrom
    alat *= BOHR_RADIUS_ANGS
    n_atoms = parse(Int, atomic_structure["nat"])

    # structure info, each column is a vector for position or lattice vector
    atom_positions = zeros(3, n_atoms)
    atom_labels = Vector{String}(undef, n_atoms)
    lattice = zeros(3, 3)

    for (i, atom) in enumerate(eachelement(findfirst("atomic_positions", atomic_structure)))
        atom_positions[:, i] = parse.(Float64, split(atom.content))
        atom_labels[i] = atom["name"]
    end
    # lattice
    for i in 1:3
        a = findfirst("cell/a$i", atomic_structure)
        lattice[:, i] = parse.(Float64, split(a.content))
    end
    # from cartesian to fractional
    atom_positions = inv(lattice) * atom_positions
    # from bohr to angstrom
    lattice *= BOHR_RADIUS_ANGS

    # reciprocal lattice
    recip_lattice = zeros(3, 3)
    for i in 1:3
        b = findfirst("basis_set/reciprocal_lattice/b$i", output)
        recip_lattice[:, i] = parse.(Float64, split(b.content))
    end
    # to 1/angstrom
    recip_lattice *= 2π / alat

    band_structure = findfirst("band_structure", output)
    n_kpts = parse(Int, findfirst("nks", band_structure).content)
    n_bands = parse(Int, findfirst("nbnd", band_structure).content)
    kpoints = zeros(3, n_kpts)
    E = zeros(n_bands, n_kpts)

    fermi_energy = parse(Float64, findfirst("fermi_energy", band_structure).content)
    # Hartree to eV
    fermi_energy *= AUTOEV

    ks_energies = findall("ks_energies", band_structure)
    for (ik, ks_energy) in enumerate(ks_energies)
        k_point = findfirst("k_point", ks_energy)
        kpoints[:, ik] = parse.(Float64, split(k_point.content))
        eigenvalues = findfirst("eigenvalues", ks_energy)
        E[:, ik] = parse.(Float64, split(eigenvalues.content))
    end
    # to 1/angstrom
    kpoints *= 2π / alat
    # from cartesian to fractional
    kpoints = inv(recip_lattice) * kpoints
    # Hartree to eV
    E *= AUTOEV

    return (; lattice, atom_positions, atom_labels, recip_lattice, kpoints, E, fermi_energy)
end
