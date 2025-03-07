
import numpy as np
import matplotlib.pyplot as plt


r_values = np.array([3.0, 3.5, 3.85, 4.0, 4.5, 5.0])  # Distances (Å)
E_total_a_u = np.array([-1054.272418, -1054.275685, -1054.275820, -1054.275773, -1054.275567, -1054.275433])  # Total energy (a.u.)
E_Ar_a_u = -527.137635  # Energy of a single argon atom (a.u.)
conversion_factor = 27.211

# Calculate binding energy in a.u.
E_binding_a_u = E_total_a_u - 2 * E_Ar_a_u

# Convert binding energy to eV
E_binding_eV = E_binding_a_u * conversion_factor

# Correcting the data and ensuring both arrays are properly aligned
r_values_corrected = np.array([3.0, 3.5, 3.85, 4.0, 4.5, 5.0])  # Explicitly redefined as 1D array
E_binding_eV_corrected = np.array([0.0776, -0.0113, -0.0150, -0.0137, -0.0081, -0.0044])  # Explicitly redefined as 1D

plt.figure(figsize=(10, 6))
plt.plot(r_values_corrected, E_binding_eV_corrected, marker='o', linestyle='-', color='purple', label="Dissociation Curve (E_binding)")
plt.axhline(0, color='red', linestyle='--', label="Zero Energy Level")

for i in range(len(r_values_corrected)):
    plt.annotate(
        f"({r_values_corrected[i]:.2f}, {E_binding_eV_corrected[i]:.4f} eV)",
        (r_values_corrected[i], E_binding_eV_corrected[i]),
        textcoords="offset points",
        xytext=(0, -20 if i % 2 == 0 else 20),
        ha='center',
        fontsize=10
    )


plt.title("Dissociation Curve for Ar$_2$ Dimer", fontsize=14)
plt.xlabel("Distance (r) [Å]", fontsize=12)
plt.ylabel("Binding Energy (E_binding) [eV]", fontsize=12)
plt.grid(alpha=0.5)
plt.legend(fontsize=10)

plt.show()

list(zip(r_values, E_binding_a_u, E_binding_eV))
manual_results = []
for i, r in enumerate(r_values):
    E_binding_manual_a_u = E_total_a_u[i] - 2 * E_Ar_a_u
    E_binding_manual_eV = E_binding_manual_a_u * conversion_factor
    manual_results.append((r, E_total_a_u[i], E_binding_manual_a_u, E_binding_manual_eV))

print("Manual Calculations for Binding Energy (E_binding):")
print(f"{'r (Å)':<10} {'E_total (a.u.)':<20} {'E_binding (a.u.)':<20} {'E_binding (eV)':<20}")
print("-" * 70)
for r, E_total, E_binding_a_u, E_binding_eV in manual_results:
    print(f"{r:<10.2f} {E_total:<20.6f} {E_binding_a_u:<20.6f} {E_binding_eV:<20.6f}")

def lennard_jones_potential(r, epsilon, sigma):
    return 4 * epsilon * ((sigma / r) ** 12 - (sigma / r) ** 6)

epsilon = 0.0103
sigma = 3.43

r_smooth = np.linspace(3.0, 5.0, 500)
lj_potential = lennard_jones_potential(r_smooth, epsilon, sigma)


plt.figure(figsize=(10, 6))
plt.plot(r_smooth, lj_potential, linestyle='--', color='blue', label="Lennard-Jones Potential")
plt.plot(r_values_corrected, E_binding_eV_corrected, marker='o', linestyle='-', color='purple', label="Dissociation Curve (E_binding)")

for i in range(len(r_values_corrected)):
    plt.annotate(
        f"({r_values_corrected[i]:.2f}, {E_binding_eV_corrected[i]:.4f} eV)",
        (r_values_corrected[i], E_binding_eV_corrected[i]),
        textcoords="offset points",
        xytext=(0, -20 if i % 2 == 0 else 20),  # Alternate above and below for readability
        ha='center',
        fontsize=10
    )

plt.axhline(-epsilon, color='red', linestyle='--', label=f"Epsilon (ε) = {-epsilon:.4f} eV")
plt.axvline(sigma, color='green', linestyle='--', label=f"Sigma (σ) = {sigma:.2f} Å")

plt.title("Lennard-Jones Potential and Dissociation Curve with Annotations", fontsize=14)
plt.xlabel("Distance (r) [Å]", fontsize=12)
plt.ylabel("Energy (E) [eV]", fontsize=12)
plt.grid(alpha=0.5)
plt.legend(fontsize=10)
plt.show()

