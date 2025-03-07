import pandas as pd
def identify_functional_groups(sdf_file_path):
    molecules = []
    current_molecule = []

    with open(sdf_file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line == "$$$$":
                molecules.append(current_molecule)
                current_molecule = []
            else:
                current_molecule.append(line)

    results = []

    for idx, molecule in enumerate(molecules):
        atom_section = []
        bond_section = []
        chemical_name = f"Molecule {idx + 1}"

        for i, line in enumerate(molecule):
            if i >= 4:
                if len(line.split()) == 16:
                    atom_section.append(line)
                elif len(line.split()) >= 3:
                    bond_section.append(line)
            if ">  <Source_ChemicalName>" in line:
                chemical_name = molecule[i + 1].strip()

        atoms = {i + 1: line.split()[3] for i, line in enumerate(atom_section)}  #Connection table oluşturarak hangi atomun hangi atomlar ile bağ yaptığını bulmaya çalıştım.
        connections = {i + 1: [] for i in range(len(atom_section))}

        for line in bond_section:
            try:
                bond_data = line.split()
                atom1 = int(bond_data[0])
                atom2 = int(bond_data[1])
                bond_type = int(bond_data[2])
                connections[atom1].append((atom2, bond_type))
                connections[atom2].append((atom1, bond_type))
            except (ValueError, IndexError):
                continue

        functional_groups = {
            "-OH": False,
            "-S=O": False,
            "-N-C=O": False,
            "-NH-C=O": False,
            "-NH2": False,
            "-O-": False,
        }

        valence = { #SDF dosyasında H atomları belirtilmediği için valans e sayısına göre H hesabı için valans hesabı yaptırdım.
            "C": 4,
            "N": 3,
            "O": 2,
            "S": 2,
        }
        for atom_id, atom_type in atoms.items():
            total_bonds = sum(bond_type for _, bond_type in connections[atom_id])
            missing_bonds = valence.get(atom_type, 0) - total_bonds #Missing bonds burada aslında H atomlarını temsil ediyor Valans değerinden yapılan bağı çıkarınca geride kalanlar H atomu olarak sayılıyor.

            if atom_type == "O":
                neighbors = connections[atom_id]
                c_bonded = any(atoms[neighbor] == "C" and bond_type == 1 for neighbor, bond_type in neighbors)
                if c_bonded and missing_bonds == 1:
                    functional_groups["-OH"] = True

            if atom_type == "N":
                neighbors = connections[atom_id]
                c_count = sum(1 for neighbor, bond_type in neighbors if atoms[neighbor] == "C" and bond_type == 1)
                if c_count == 1 and missing_bonds == 2:
                    functional_groups["-NH2"] = True

                if c_count >= 1:
                    for neighbor, bond_type in neighbors:
                        if atoms[neighbor] == "C":
                            carbon_neighbors = connections[neighbor]
                            if any(atoms[c_neighbor] == "O" and c_bond_type == 2 for c_neighbor, c_bond_type in carbon_neighbors):
                                if missing_bonds == 1:
                                    functional_groups["-NH-C=O"] = True
                                elif missing_bonds == 0:
                                    functional_groups["-N-C=O"] = True

            if atom_type == "S":
                neighbors = connections[atom_id]
                o_double_bonded = sum(1 for neighbor, bond_type in neighbors if atoms[neighbor] == "O" and bond_type == 2)
                if o_double_bonded >= 1:
                    functional_groups["-S=O"] = True

            if atom_type == "O":
                neighbors = connections[atom_id]
                c_count = sum(1 for neighbor, bond_type in neighbors if atoms[neighbor] == "C" and bond_type == 1)
                if c_count == 2:
                    functional_groups["-O-"] = True

        results.append({
            "Chemical Name": chemical_name,
            **functional_groups
        })

    return results

sdf_file_path = "sample.sdf"
results = identify_functional_groups(sdf_file_path)
df = pd.DataFrame(results)
output_csv = "functional_group_analysis"
df.to_csv(output_csv, index=False)
print(df)
