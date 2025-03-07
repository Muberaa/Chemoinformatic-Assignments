import numpy as np
import matplotlib.pyplot as plt

# Define the Lennard-Jones potential function
def lennard_jones_potential(r, epsilon, sigma):
    return 4 * epsilon * ((sigma / r) ** 12 - (sigma / r) ** 6)

# Given parameters
epsilon = 0.0103  # Depth of the potential well in eV
sigma = 3.43  # Distance where the potential is at its minimum (Å)

# Define r values
r_given = np.array([3.0, 3.5, 3.85, 4.0, 4.5, 5.0])  # Specific r values for calculation
r_values_lj = np.linspace(2.5, 6.0, 500)  # Smooth range for Lennard-Jones curve

# Calculate Lennard-Jones potential
lj_potential_given = lennard_jones_potential(r_given, epsilon, sigma)  # For given r values
lj_potential_full = lennard_jones_potential(r_values_lj, epsilon, sigma)  # For smooth curve

# Manual calculation for the given r values
manual_calculations = []
for r in r_given:
    value = 4 * epsilon * ((sigma / r) ** 12 - (sigma / r) ** 6)
    manual_calculations.append(value)

# Print manual calculations
print("Manual Calculations:")
for r, e in zip(r_given, manual_calculations):
    print(f"r = {r:.2f} Å, E(r) = {e:.4f} eV")

# Plot Lennard-Jones potential for smooth and given points
plt.figure(figsize=(12, 8))

# Smooth Lennard-Jones curve
plt.plot(r_values_lj, lj_potential_full, linestyle='-', color='blue', label="Lennard-Jones Potential (Smooth)")

# Points for given r values
plt.plot(r_given, lj_potential_given, marker='o', linestyle='-', color='orange', label="Given Points")


for i in range(len(r_given)):
    plt.annotate(
        f"({r_given[i]:.2f}, {lj_potential_given[i]:.4f})",
        (r_given[i], lj_potential_given[i]),
        textcoords="offset points",
        xytext=(0, -15 if i % 2 == 0 else 15),  # Alternate above and below for clarity
        ha='center',
        fontsize=10
    )


plt.axhline(-epsilon, color='red', linestyle='--', label=f"Epsilon (ε) = {-epsilon:.4f} eV")
plt.axvline(sigma, color='green', linestyle='--', label=f"Sigma (σ) = {sigma:.2f} Å")
plt.title("Lennard-Jones Potential for Given and Smooth r Values", fontsize=16)
plt.xlabel("Distance (r) [Å]", fontsize=14)
plt.ylabel("Potential Energy (E) [eV]", fontsize=14)
plt.grid(alpha=0.5)
plt.legend(fontsize=12)

plt.show()