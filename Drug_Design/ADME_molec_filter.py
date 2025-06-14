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
