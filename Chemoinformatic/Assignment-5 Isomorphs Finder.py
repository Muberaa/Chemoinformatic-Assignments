import itertools
import numpy as np
def parse_mol_file(file_path):
    atoms = []
    bonds = []

    with open(file_path, 'r') as file:
        lines = file.readlines()

    atom_count = int(lines[3][:3].strip())
    bond_count = int(lines[3][3:6].strip())

    for i in range(4, 4 + atom_count):
        parts = lines[i].split()
        atoms.append(parts[3])

    for i in range(4 + atom_count, 4 + atom_count + bond_count):
        parts = lines[i].split()
        bonds.append((int(parts[0]), int(parts[1]), int(parts[2])))
    return atoms, bonds
def build_adjacency_matrix(atoms, bonds): #Atom sayısına göre bağ matrisi oluşturup atomlar için bağ tipi belirledim.
    n = len(atoms)
    adj_matrix = np.zeros((n + 1, n + 1), dtype=int)

    for bond in bonds:
        i, j, bond_type = bond
        adj_matrix[i][j] = bond_type
        adj_matrix[j][i] = bond_type

    return adj_matrix[1:, 1:]
def find_total_isomorphs(atoms, bonds):
    adj_matrix = build_adjacency_matrix(atoms, bonds)
    n = len(atoms)

    total_isomorphs = len(list(itertools.permutations(range(n)))) #Atom sayılarına göre bir permütasyon oluşturup N! formülünden toplam permütasyon sayısını bulmak istedim.
    return total_isomorphs
def main(file_path):
    atoms, bonds = parse_mol_file(file_path)
    print(f"Parsed Atoms: {atoms}")
    print(f"Parsed Bonds: {bonds}")

    total_isomorphs = find_total_isomorphs(atoms, bonds)
    print(f"Total Isomorphisms: {total_isomorphs}")

if __name__ == "__main__":
    file_path = "acetic_acid.mol.txt"
    try:
        main(file_path)
    except Exception as e:
        print(f"Error: {e}")