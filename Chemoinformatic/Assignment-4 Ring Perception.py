import pandas as pd

def parse_mol_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    atom_count = int(lines[3][:3].strip())
    bond_count = int(lines[3][3:6].strip())

    atoms = []
    for i in range(4, 4 + atom_count):
        line = lines[i].strip()
        atoms.append({"index": i - 4, "type": line[31:34].strip()})

    bonds = []
    for i in range(4 + atom_count, 4 + atom_count + bond_count):
        line = lines[i].strip().split()
        atom1 = int(line[0]) - 1
        atom2 = int(line[1]) - 1
        bond_type = int(line[2])
        bonds.append({"atom1": atom1, "atom2": atom2, "type": bond_type})

    return atoms, bonds

def build_graph(atoms, bonds):
    graph = {atom["index"]: [] for atom in atoms} #
    for bond in bonds:
        graph[bond["atom1"]].append((bond["atom2"], bond["type"]))
        graph[bond["atom2"]].append((bond["atom1"], bond["type"]))
    return graph

def find_cycles(graph): # DFS ile atomların komşuluğuna göre bir döngü oluşturup bu döngüleri benzersizleştirerek halka yapı olup olmadığını kontrol ettim
    def dfs(node, visited, stack, parent):
        visited[node] = True
        stack.append(node)

        for neighbor, _ in graph[node]:
            if not visited.get(neighbor, False):
                if dfs(neighbor, visited, stack, node):
                    return True
            elif neighbor in stack and neighbor != parent:
                cycle_start_idx = stack.index(neighbor)
                cycles.append(stack[cycle_start_idx:])

        stack.pop()
        return False

    visited = {}
    stack = []
    cycles = []

    for node in graph:
        if not visited.get(node, False):
            dfs(node, visited, stack, -1)

    unique_cycles = []
    for cycle in cycles:
        if sorted(cycle) not in unique_cycles:
            unique_cycles.append(sorted(cycle))
    return unique_cycles

def is_aromatic_ring(atoms, bonds, ring): #Aromatik olup olmadıklarını belirlemek için halka için bond type ve atom sayısının kurallarını girdim
    if len(ring) != 6:
        return False

    bond_map = {(bond["atom1"], bond["atom2"]): bond["type"] for bond in bonds}
    bond_map.update({(bond["atom2"], bond["atom1"]): bond["type"] for bond in bonds})

    bond_types = []
    for i in range(len(ring)):
        atom1 = ring[i]
        atom2 = ring[(i + 1) % len(ring)]
        bond_types.append(bond_map.get((atom1, atom2), 0))

    alternating_pattern = [1, 2] * 3
    return bond_types == alternating_pattern or bond_types == alternating_pattern[::-1]


def analyze_molecule(file_path):
    atoms, bonds = parse_mol_file(file_path)
    molecule_graph = build_graph(atoms, bonds)

    cycles = find_cycles(molecule_graph)

    aromatic_rings = []
    ring_atoms = []

    for cycle in cycles:
        ring_atoms.append(cycle)
        if is_aromatic_ring(atoms, bonds, cycle):
            aromatic_rings.append(cycle)

    results = {
        "Has Ring": len(cycles) > 0,
        "Total Rings": len(cycles),
        "Ring Atoms": ring_atoms,
        "Aromatic Rings": aromatic_rings
    }

    return results


file_path = "phenylalanine.mol.txt"
results = analyze_molecule(file_path)

results_df = pd.DataFrame({
    "Description": ["Has Ring", "Total Rings", "Aromatic Rings"],
    "Value": [results["Has Ring"], results["Total Rings"], len(results["Aromatic Rings"])]
})

print(results_df)
if results["Ring Atoms"]:
    print("List of Aromatic Ring Atoms:", results["Aromatic Rings"])

