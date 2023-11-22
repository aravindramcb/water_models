import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
# Data for the properties and their corresponding values
properties = ["Density", "Coordination \nnumber", "Surface\n tension",
              "Dielectric \nconstant", "Self-diffusion\n coefficient",
              "Solvation \nfree \nenergies \nof methane", "Dipole \nmoment"]

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
sns.set_theme(style='white', context="paper", font_scale=1.3)
fig, (ax1,ax2) = plt.subplots(2,1, figsize=(10,5), dpi=300, sharex=True,
                              gridspec_kw={'height_ratios': [1.5, 5]})
ax1.set_position([0.1, 0.3, 0.8, 0.2])  # [left, bottom, width, height]
fig.subplots_adjust(hspace=0.15)

# Top part
ax1.bar(x - width, relative_errors["OPC"], width, label='OPC')
ax1.bar(x, relative_errors["TIP3P"], width, label='TIP3P')
ax1.bar(x + width, relative_errors["TIP4P_Ew"], width, label='TIP4P-Ew')

# Bottom part
ax2.bar(x - width, relative_errors["OPC"], width, label='OPC')
ax2.bar(x, relative_errors["TIP3P"], width, label='TIP3P')
ax2.bar(x + width, relative_errors["TIP4P_Ew"], width, label='TIP4P-Ew')

# Zoom in
ax1.set_ylim(100,160)
ax2.set_ylim(0,40)

# hide spines between axes
ax1.spines.bottom.set_visible(False)
ax2.spines.top.set_visible(False)

ax1.tick_params(labeltop=False)  # don't put tick labels at the top
ax2.xaxis.tick_bottom()
ax1.xaxis.tick_top()
# cutout lines
d = .5  # proportion of vertical to horizontal extent of the slanted line
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
              linestyle="none", color='k', mec='k', mew=1, clip_on=False)
ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)

ax2.set_ylabel('Relative Error (%)')
ax2.set_xticks(x)
ax2.set_xticklabels(properties, rotation=0)
ax1.legend(loc='upper left',ncol=3)

plt.tight_layout()
# plt.show()

plt.savefig("/home/aravind/PhD_local/dean/figures/main_images/water_properties.png")
