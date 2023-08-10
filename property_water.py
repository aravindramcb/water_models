import matplotlib.pyplot as plt
import numpy as np

# Data for the properties and their corresponding values
properties = ["Density (g/cm3)", "Coordination number", "Surface tension (mN/m)", 
              "Dielectric constant", "Self-diffusion coefficient (× 10-5 cm2 /s)", 
              "Δgmethane (kJ/mol)", "Dipole moment"]

# Values for the different water models
OPC = [0.997, 5.1971, 70.1, 78.00, 2.27, 9.34, 2.48]
TIP3P = [0.98, 6.239, 47, 95.00, 5.72, 8.51, 2.35]
TIP4P_Ew = [0.996, 4.69, 59.2, 65.00, 2.54, 8.8, 2.32]
Experiment = [0.997, 4.7, 71.99, 78.30, 2.30, 8.803, 1.84]

# Calculating the relative errors for each property
def calculate_relative_error(experimental, theoretical):
    return abs((experimental - theoretical) / experimental) * 100

relative_errors = {
    "OPC": [calculate_relative_error(Experiment[i], OPC[i]) for i in range(len(properties))],
    "TIP3P": [calculate_relative_error(Experiment[i], TIP3P[i]) for i in range(len(properties))],
    "TIP4P_Ew": [calculate_relative_error(Experiment[i], TIP4P_Ew[i]) for i in range(len(properties))]
}

# Creating the clustered column chart
x = np.arange(len(properties))  # X-axis values for properties
width = 0.2  # Width of each column

fig, ax = plt.subplots(figsize=(8,4),dpi=300)
rects1 = ax.bar(x - width, relative_errors["OPC"], width, label='OPC')
rects2 = ax.bar(x, relative_errors["TIP3P"], width, label='TIP3P')
rects3 = ax.bar(x + width, relative_errors["TIP4P_Ew"], width, label='TIP4P-Ew')

ax.set_xlabel('Properties')
ax.set_ylabel('Relative Error (%)')
ax.set_title('Relative Error of Water Models')
ax.set_xticks(x)
ax.set_xticklabels(properties, rotation=45)
ax.legend()

fig.tight_layout()
plt.savefig("/home/aravind/PhD_local/dean/figures/main_images/water_properties.png")
