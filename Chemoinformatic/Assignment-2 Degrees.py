
with open("acetic_acid.mol.txt", "r") as file:
    acetic_acid = file.read()

with open("phenylalanine.mol.txt", "r") as file:
    phenylalanine = file.read()

def parse_mol_file(mol_content):

    lines = mol_content.splitlines()


    atom_count = int(lines[3][:3].strip())
    bond_count = int(lines[3][3:6].strip())


    bond_block = lines[4 + atom_count:4 + atom_count + bond_count]


    degrees = [0] * atom_count


    for bond in bond_block:
        atom1 = int(bond[:3].strip()) - 1
        atom2 = int(bond[3:6].strip()) - 1
        degrees[atom1] += 1
        degrees[atom2] += 1

    return degrees


acetic_acid_degrees = parse_mol_file(acetic_acid)
phenylalanine_degrees = parse_mol_file(phenylalanine)


print("Acetic Acid Atom Degrees:", acetic_acid_degrees)
print("Phenylalanine Atom Degrees:", phenylalanine_degrees)