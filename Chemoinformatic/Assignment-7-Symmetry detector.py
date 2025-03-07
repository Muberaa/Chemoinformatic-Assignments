import numpy as np
from scipy.spatial.distance import cdist
import itertools

def parse_xyz(file_content):
    molecules = []
    lines = file_content.strip().split("\n")

    idx = 0
    while idx < len(lines):
        num_atoms = int(lines[idx].strip())
        header = lines[idx + 1].strip()
        molecule_atoms = []
        for i in range(idx + 2, idx + 2 + num_atoms):
            parts = lines[i].split()
            element = parts[0]
            coords = np.array([float(x) for x in parts[1:]])
            molecule_atoms.append((element, coords))
        molecules.append((header, molecule_atoms))
        idx += 2 + num_atoms

    return molecules

def are_symmetric(molecule1, molecule2):
    coords1 = np.array([atom[1] for atom in molecule1])
    coords2 = np.array([atom[1] for atom in molecule2])


    dist_matrix1 = cdist(coords1, coords1)
    dist_matrix2 = cdist(coords2, coords2)

    return np.allclose(np.sort(dist_matrix1.flatten()), np.sort(dist_matrix2.flatten()))

def find_symmetric_molecules(molecules):
    symmetric_map = {i + 1: set() for i in range(len(molecules))}
    for i, j in itertools.combinations(range(len(molecules)), 2):
        header1, atoms1 = molecules[i]
        header2, atoms2 = molecules[j]
        if are_symmetric(atoms1, atoms2):
            symmetric_map[i + 1].add(j + 1)
            symmetric_map[j + 1].add(i + 1)
    return symmetric_map


file_path = "Pyr_Pyr_chemo.xyz"
with open(file_path, "r") as file:
    file_content = file.read()

molecules = parse_xyz(file_content)
symmetric_map = find_symmetric_molecules(molecules)


print("Symmetric Molecule Relationships:")
for molecule, symmetric_molecules in symmetric_map.items():
    if symmetric_molecules:
        symmetric_list = ", ".join(map(str, sorted(symmetric_molecules)))
        print(f"Molecule {molecule} is symmetric with Molecules: {symmetric_list}")