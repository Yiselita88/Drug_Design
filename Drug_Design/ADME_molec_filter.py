from pathlib import Path
import math

import numpy as np
import pandas as pd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw, PandasTools
from matplotlib.lines import Line2D

import openpyxl
from openpyxl import Workbook
from openpyxl.drawing.image import Image as XLImage
from io import BytesIO

HERE = Path(__file__).resolve()
# print(HERE)
DATA = HERE/"data"

# Define & visualize molecules
smiles = [
    "CCC1C(=O)N(CC(=O)N(C(C(=O)NC(C(=O)N(C(C(=O)NC(C(=O)NC(C(=O)N(C(C(=O)N(C(C(=O)N(C(C(=O)N(C(C(=O)N1)C(C(C)CC=CC)O)C)C(C)C)C)CC(C)C)C)CC(C)C)C)C)C)CC(C)C)C)C(C)C)CC(C)C)C)C",
    "CN1CCN(CC1)C2=C3C=CC=CC3=NC4=C(N2)C=C(C=C4)C",
    "CC1=C(C(CCC1)(C)C)C=CC(=CC=CC(=CC=CC=C(C)C=CC=C(C)C=CC2=C(CCCC2(C)C)C)C)C",
    "CCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)O",
]
names = ["cyclosporine", "clozapine", "beta-carotene", "cannabidiol"]
molecules = pd.DataFrame({"name": names, "smiles": smiles})
PandasTools.AddMoleculeColumnToFrame(molecules, "smiles")
# print(molecules)

# # Optional: for saving the dataframe with the corresponding figures:
# molecules["Mol"] = molecules["smiles"].apply(Chem.MolFromSmiles)
# wb = Workbook()
# ws = wb.active
# ws.append(["Name", "SMILES", "Sructure"])

# # add rows with images
# for i, row in molecules.iterrows():
#     mol_img =  Draw.MolToImage(row["Mol"], size=(150, 150))
#     img_io = BytesIO()
#     mol_img.save(img_io, format='PNG')
#     img_io.seek(0)

#     img  = XLImage(img_io)
#     img.width = 100
#     img.height = 100

#     ws.append([row["name"], row["smiles"], ""])
#     ws.add_image(img, f"c{i+2}")

# wb.save("molecules_with_images.xlsx")

molecules["molecular_weight"] = molecules["ROMol"].apply(Descriptors.ExactMolWt)
molecules["n_hba"] = molecules["ROMol"].apply(Descriptors.NumHAcceptors)
molecules["n_hbd"] = molecules["ROMol"].apply(Descriptors.NumHDonors)
molecules["logp"] = molecules["ROMol"].apply(Descriptors.MolLogP)
molecules["color"] = ["red", "green", "blue", "cyan"]
print(molecules[["molecular_weight", "n_hba", "n_hbd", "logp"]])

ro5_properties = {
    "molecular_weight": (500, "molecular weight (Da)"),
    "n_hba": (10, "# HBA"),
    "n_hbd": (5, "# HBD"),
    "logp": (5, "logP"),
}

# Start 1x4 plot frame
fig, axes = plt.subplots(figsize=(10, 2.5), nrows=1, ncols=4)
x = np.arange(1, len(molecules) + 1)
colors = ["red", "green", "blue", "cyan"]

# Create subplots
for index, (key, (threshold, title)) in enumerate(ro5_properties.items()):
    axes[index].bar([1, 2, 3, 4], molecules[key], color=colors)
    axes[index].axhline(y=threshold, color="black", linestyle="dashed")
    axes[index].set_title(title)
    axes[index].set_xticks([])

# Add legend
legend_elements = [
    mpatches.Patch(color=row["color"], label=row["name"]) for index, row in molecules.iterrows()
]
legend_elements.append(Line2D([0], [0], color="black", ls="dashed", label="Threshold"))
fig.legend(handles=legend_elements, bbox_to_anchor=(1.2, 0.8))

# Fit subplots and legend into figure
# plt.tight_layout()
# plt.savefig('Ro5_molecules.png')

# Investigate compliance with Ro5
def calculate_Ro5_properties(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    molecular_weight = Descriptors.ExactMolWt(molecule)
    n_hba = Descriptors.NumHAcceptors(molecule)
    n_hbd = Descriptors.NumHDonors(molecule)
    logp = Descriptors.MolLogP(molecule)

    conditions = [molecular_weight <= 500, n_hba <=10, n_hbd <=5, logp <= 5]
    ro5_fulfilled = sum(conditions) >= 3 
    return pd.Series([molecular_weight, n_hba, n_hbd, logp, ro5_fulfilled], index=["molecular_weight", "n_hba", "n_hbd", "logp", "ro5_fulfilled"],)

for name, smiles in zip(molecules["name"], molecules["smiles"]):
    print(f"Ro5 fulfilled for {name}: {calculate_Ro5_properties(smiles)['ro5_fulfilled']}")

# # Apply Ro5 ti EGFR dataset
molecules = pd.read_csv("EGFR_compounds.csv", index_col=0)
print(molecules.shape)
print(molecules.head())
ro5_properties = molecules["smiles"].apply(calculate_Ro5_properties)
print(ro5_properties.head())

# concatenate molecules with Ro5 data
molecules = pd.concat([molecules, ro5_properties], axis=1)
print(molecules.head())

molecules_ro5_fulfilled = molecules[molecules["ro5_fulfilled"]]
molecules_ro5_violated = molecules[~molecules["ro5_fulfilled"]]

print(f"# compounds in unfiltered data set: {molecules.shape[0]}")
print(f"# compounds in filtered data set: {molecules_ro5_fulfilled.shape[0]}")
print(f"# compounds not compliant with the Ro5: {molecules_ro5_violated.shape[0]}")

molecules_ro5_fulfilled.to_csv("EGFR_compounds_lipinski.csv")
molecules_ro5_fulfilled.head()

# # Calculate statistics on Ro5 properties
def calculate_mean_std(dataframe):
    stats = dataframe.describe()
    stats = stats.T
    stats = stats[["mean", "std"]]
    return stats

molecules_ro5_fulfilled_stats = calculate_mean_std(molecules_ro5_fulfilled[["molecular_weight", "n_hba", "n_hbd", "logp"]])
print(molecules_ro5_fulfilled_stats)

molecules_ro5_violated_stats = calculate_mean_std(molecules_ro5_violated[["molecular_weight", "n_hba", "n_hbd", "logp"]])
print(molecules_ro5_violated_stats)

# # helper function for better scaling and further radar plot
def _scale_by_thresholds(stats, thresholds, scaled_treshold):
    for property_name in stats.index:
        if property_name not in thresholds.keys():
            raise KeyError(f"Add property '{property_name}' to scaling variable")
    
    stats_scaled = stats.apply(lambda x: x / thresholds[x.name] * scaled_treshold, axis=1)
    return stats_scaled

def _define_radial_axes_angles(n_axes):
    x_angles = [i / float(n_axes) * 2 * math.pi for i in range(n_axes)]
    x_angles += x_angles[:1]
    return x_angles

def plot_radar(
    y,
    thresholds,
    scaled_threshold,
    properties_labels,
    y_max=None,
    output_path=None,
):
    """
    Plot a radar chart based on the mean and standard deviation of a data set's properties.

    Parameters
    ----------
    y : pd.DataFrame
        Dataframe with "mean" and "std" (columns) for each physicochemical property (rows).
    thresholds : dict of str: int
        Thresholds defined for each property.
    scaled_threshold : int or float
        Scaled thresholds across all properties.
    properties_labels : list of str
        List of property names to be used as labels in the plot.
    y_max : None or int or float
        Set maximum y value. If None, let matplotlib decide.
    output_path : None or pathlib.Path
        If not None, save plot to file.
    """

    # Define radial x-axes angles -- uses our helper function!
    x = _define_radial_axes_angles(len(y))
    # Scale y-axis values with respect to a defined threshold -- uses our helper function!
    y = _scale_by_thresholds(y, thresholds, scaled_threshold)
    # Since our chart will be circular we append the first value of each property to the end
    y = pd.concat([y, y.head(1)])

    # Set figure and subplot axis
    plt.figure(figsize=(6, 6))
    ax = plt.subplot(111, polar=True)

    # Plot data
    ax.fill(x, [scaled_threshold] * len(x), "cornflowerblue", alpha=0.2)
    ax.plot(x, y["mean"], "b", lw=3, ls="-")
    ax.plot(x, y["mean"] + y["std"], "orange", lw=2, ls="--")
    ax.plot(x, y["mean"] - y["std"], "orange", lw=2, ls="-.")

    # From here on, we only do plot cosmetics
    # Set 0° to 12 o'clock
    ax.set_theta_offset(math.pi / 2)
    # Set clockwise rotation
    ax.set_theta_direction(-1)

    # Set y-labels next to 180° radius axis
    ax.set_rlabel_position(180)
    # Set number of radial axes' ticks and remove labels
    plt.xticks(x, [])
    # Get maximal y-ticks value
    if not y_max:
        y_max = int(ax.get_yticks()[-1])
    # Set axes limits
    plt.ylim(0, y_max)
    # Set number and labels of y axis ticks
    plt.yticks(
        range(1, y_max),
        ["5" if i == scaled_threshold else "" for i in range(1, y_max)],
        fontsize=16,
    )

    # Draw ytick labels to make sure they fit properly
    # Note that we use [:1] to exclude the last element which equals the first element (not needed here)
    for i, (angle, label) in enumerate(zip(x[:-1], properties_labels)):
        if angle == 0:
            ha = "center"
        elif 0 < angle < math.pi:
            ha = "left"
        elif angle == math.pi:
            ha = "center"
        else:
            ha = "right"
        ax.text(
            x=angle,
            y=y_max + 1,
            s=label,
            size=16,
            horizontalalignment=ha,
            verticalalignment="center",
        )

    # Add legend relative to top-left plot
    labels = ("mean", "mean + std", "mean - std", "rule of five area")
    ax.legend(labels, loc=(1.1, 0.7), labelspacing=0.3, fontsize=16)

    # Save plot - use bbox_inches to include text boxes
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches="tight", transparent=True)

    plt.savefig("radar_plot.png")

thresholds = {"molecular_weight": 500, "n_hba": 10, "n_hbd": 5, "logp": 5}
scaled_threshold = 5
properties_labels = [
    "Molecular weight (Da) / 100",
    "# HBA / 2",
    "# HBD",
    "LogP",
]
y_max = 8

plot_radar(molecules_ro5_fulfilled_stats, thresholds, scaled_threshold, properties_labels, y_max)