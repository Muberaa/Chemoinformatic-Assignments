
class MolecularGraph:
    def __init__(self):
        self.atoms = []
        self.bonds = {}
        self.atom_symbols = {}
    def add_atom(self, atom_id, symbol="C"):
        if atom_id not in self.bonds:
            self.bonds[atom_id] = []
            self.atoms.append(atom_id)
            self.atom_symbols[atom_id] = symbol

    def add_bond(self, atom1, atom2):
        self.bonds[atom1].append(atom2)
        self.bonds[atom2].append(atom1)

    def degree(self, atom_id):
        return len(self.bonds[atom_id])

def parse_mol_file(filepath):
    graph = MolecularGraph()
    with open(filepath, 'r') as file:
        lines = file.readlines()

    atom_count, bond_count = map(int, lines[3].split()[:2])

    for i in range(4, 4 + atom_count):
        parts = lines[i].split()
        atom_symbol = parts[3]
        graph.add_atom(i - 3, atom_symbol)

    for i in range(4 + atom_count, 4 + atom_count + bond_count):
        atom1, atom2, _ = map(int, lines[i].split()[:3])
        graph.add_bond(atom1, atom2)
    return graph

def morgan_algorithm(graph, max_iterations=10): #burada atom derecelerine göre bir invariant atadım
   #bunu da bond table ile oluşturmuştum zaten daha sonra oluşan invariantı yine komşu atomların invariantlarını toplayıp sayı değerinin artmadığı yere kadar döngü oluşturdum
   #stabil olmasından emin olabilmek için iterasyon sayısını max olarak tutturdum.

    invariants = {atom: graph.degree(atom) for atom in graph.atoms}
    history = {atom: [graph.degree(atom)] for atom in graph.atoms}

    for _ in range(max_iterations):
        new_invariants = {}
        for atom in graph.atoms:
            new_value = sum(invariants[neighbor] for neighbor in graph.bonds[atom])
            new_invariants[atom] = new_value
            history[atom].append(new_value)

        if new_invariants == invariants:
            break
        invariants = new_invariants

    sorted_atoms = sorted(invariants.items(), key=lambda x: x[1])
    ranks = {}
    current_rank = 1
    last_value = None
    for atom, value in sorted_atoms:
        if value != last_value:
            ranks[atom] = current_rank
            last_value = value
            current_rank += 1
        else:
            ranks[atom] = current_rank - 1

    unique_counts = []
    for i in range(len(next(iter(history.values())))):
        unique_values = {history[atom][i] for atom in history}
        unique_counts.append(len(unique_values))

    return ranks, history, unique_counts

def format_output(title, ranks, history, unique_counts, graph):
    output = [f"{title}:"]
    output.append("Atom ID | Symbol | Rank | Invariant History")
    output.append("--------|--------|------|------------------")
    for atom in sorted(graph.atoms):
        symbol = graph.atom_symbols[atom]
        rank = ranks[atom]
        history_str = " -> ".join(map(str, history[atom]))
        output.append(f"{atom:<8}| {symbol:<6}| {rank:<5}| {history_str}")

    output.append("\nUnique Invariant Counts per Iteration:")
    output.append(" -> ".join(map(str, unique_counts)))
    return "\n".join(output)

phenylalanine_graph = parse_mol_file('phenylalanine.mol.txt')
acetic_acid_graph = parse_mol_file('acetic_acid.mol.txt')

phenylalanine_ranks, phenylalanine_history, phenylalanine_unique_counts = morgan_algorithm(phenylalanine_graph)
acetic_acid_ranks, acetic_acid_history, acetic_acid_unique_counts = morgan_algorithm(acetic_acid_graph)

phenylalanine_output = format_output("Phenylalanine Atom Ranks", phenylalanine_ranks, phenylalanine_history,
                                     phenylalanine_unique_counts, phenylalanine_graph)
acetic_acid_output = format_output("Acetic Acid Atom Ranks", acetic_acid_ranks, acetic_acid_history,
                                   acetic_acid_unique_counts, acetic_acid_graph)

print(phenylalanine_output)
print("\n" + acetic_acid_output)