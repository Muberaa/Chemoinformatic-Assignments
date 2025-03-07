sdf_file_path = "sample.sdf"

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

for molecule in molecules:
    atom_count = 0
    molecular_weight = None
    chemical_name = "Unknown"


    if len(molecule) > 3:
        counts_line = molecule[3]
        parts = counts_line.split()
        atom_count = int(parts[0]) if parts else 0


    for i, line in enumerate(molecule):

        if "STRUCTURE_MolecularWeight" in line and molecular_weight is None:
            try:
                molecular_weight = float(molecule[i + 1].strip())
            except ValueError:
                molecular_weight = None


        elif "Source_ChemicalName" in line and chemical_name == "Unknown":
            chemical_name = molecule[i + 1].strip()


    results.append({
        "Total Atom Number": atom_count,
        "Molecular Weight": molecular_weight,
        "Source Chemical Name": chemical_name
    })


print(f"Total molecules in file: {len(results)}\n")


for i, data in enumerate(results, start=1):
    print(f"Molecule {i}:")
    print(f"Total Atom Number: {data['Total Atom Number']}")
    print(f"Molecular Weight: {data['Molecular Weight']:.2f}" if data['Molecular Weight']
          else "Molecular Weight: Not available")
    print(f"Source Chemical Name: {data['Source Chemical Name']}\n")



