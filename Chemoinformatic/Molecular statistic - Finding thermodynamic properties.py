import numpy as np
import pandas as pd

def calculate_properties_from_table(data, temperatures):
    # Constants
    h = 6.62607015e-34
    k = 1.380649e-23
    R = 8.314
    P = 1e5
    NA = 6.02214076e23
    c = 2.99792458e10

    def translational_entropy(T, m):
        q_trans = ((2 * np.pi * m * k * T) / h**2)**(3 / 2) * (k * T / P)
        return np.log(q_trans) * R

    def electronic_partition_function(T, energy_levels):
        q_elec = sum(g * np.exp(-((E * c * h) / (k * T))) for E, g in energy_levels)
        return q_elec

    results = []

    for atom in data:
        name = atom["name"]
        m = atom["mass"] / 1000 / NA
        energy_levels = atom.get("energy_levels", [(0, 1)])

        atom_results = {"name": name, "properties": []}

        for T in temperatures:
            S_trans = translational_entropy(T, m)


            q_elec = electronic_partition_function(T, energy_levels)

            E_trans = (3 / 2) * R * T
            E_elec = sum((g * E * c * h * np.exp(-((E * c * h) / (k * T)))) for E, g in energy_levels) / q_elec
            E = E_trans + E_elec  # Total internal energy

            # Enthalpy (H)
            H = E + R * T

            # Entropy (S)
            S_elec = R * (np.log(q_elec))
            S = S_trans + S_elec + R * 5 / 2  # Total entropy

            # Gibbs free energy (G)
            G = H - T * S

            # Heat capacity
            Cv = (3 / 2) * R  # Translational contribution only

            atom_results["properties"].append({
                "T": T,
                "E (kJ/mol)": E / 1000,  # Convert to kJ/mol
                "H (kJ/mol)": H / 1000,  # Convert to kJ/mol
                "S (J/(mol*K))": S,  # J/(mol*K)
                "G (kJ/mol)": G / 1000,  # Convert to kJ/mol
                "Cv (J/(mol*K))": Cv  # J/(mol*K)
            })

        results.append(atom_results)

    return results

# Data for atoms from the assignment table
# Includes energy levels (E0, E1, g_e) explicitly
data = [
    {"name": "H", "mass": 1.008, "energy_levels": [(0, 2), (82259, 2)]},
    {"name": "He", "mass": 4.0026, "energy_levels": [(0, 1), (159856, 3)]},
    {"name": "Li", "mass": 6.94, "energy_levels": [(0, 2), (14904, 2)]},
    {"name": "Be", "mass": 9.0122, "energy_levels": [(0, 1), (21978, 1)]},
    {"name": "B", "mass": 10.81, "energy_levels": [(0, 2), (15.3, 4)]},
    {"name": "C", "mass": 12.011, "energy_levels": [(0, 1), (16.4, 3), (43.4, 5)]},
    {"name": "N", "mass": 14.007, "energy_levels": [(0, 4), (19225, 6)]},
    {"name": "O", "mass": 15.999, "energy_levels": [(0, 5), (158.3, 3), (226.5, 1)]},
    {"name": "F", "mass": 18.998, "energy_levels": [(0, 4), (404, 2)]},
    {"name": "Ne", "mass": 20.180, "energy_levels": [(0, 1), (134042, 5)]},
    {"name": "Na", "mass": 22.990, "energy_levels": [(0, 2), (16956, 2)]},
    {"name": "Ar", "mass": 39.948, "energy_levels": [(0, 1), (93144, 5)]}
]

temperatures = [300, 500, 1000, 2000, 3000, 5000]

results = calculate_properties_from_table(data, temperatures)

for atom in results:
    atom_name = atom["name"]
    rows = []
    for prop in atom["properties"]:
        rows.append({
            "T (K)": prop["T"],
            "molar E (kJ/mol)": prop["E (kJ/mol)"],
            "molar H (kJ/mol)": prop["H (kJ/mol)"],
            "molar S (J/(mol*K))": prop["S (J/(mol*K))"],
            "molar G (kJ/mol)": prop["G (kJ/mol)"],
            "molar Cv (J/(mol*K))": prop["Cv (J/(mol*K))"]
        })
    df_atom = pd.DataFrame(rows)
    output_file = f"thermodynamic_properties_{atom_name}.csv"
    df_atom.to_csv(output_file, index=False)
    print(f"Results for {atom_name} have been saved to '{output_file}'")