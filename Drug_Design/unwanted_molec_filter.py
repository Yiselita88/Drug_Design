from pathlib import Path

import pandas as pd
import tqdm.auto as tqdm 
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

egfr_data = pd.read_csv("EGFR_compounds_lipinski.csv", index_col=0)
print("Dataframe shape:", egfr_data.shape)
egfr_data.drop(columns=["molecular_weight", "n_hbd", "n_hba", "logp"], inplace=True)
# print(egfr_data.head())

PandasTools.AddMoleculeColumnToFrame(egfr_data, smilesCol="smiles")
Chem.Draw.MolsToGridImage(list(egfr_data.head(3).ROMol), legends=list(egfr_data.head(3).molecule_chembl_id))
# img = Chem.Draw.MolsToGridImage(list(egfr_data.head(3).ROMol), legends=list(egfr_data.head(3).molecule_chembl_id))
# img.save("data_T003.png")

# # Filter for PAINS
substructures = pd.read_csv("unwanted_substructures.csv", sep="\s+")
substructures["rdkit_molecule"] = substructures.smarts.apply(Chem.MolFromSmarts)
print("Number of unwanted substructures in collection:", len(substructures))

Chem.Draw.MolsToGridImage(mols=substructures.rdkit_molecule.tolist()[2:5], legends=substructures.name.tolist()[2:5])
# img2 = Chem.Draw.MolsToGridImage(mols=substructures.rdkit_molecule.tolist()[2:5], legends=substructures.name.tolist()[2:5])
# img2.save("substructures.png")

# # Search the filtered DF for matches of unwanted substructures
matches = []
clean = []
for index, row in tqdm.tqdm(egfr_data.iterrows(), total=egfr_data.shape[0]):
    molecule = Chem.MolFromSmiles(row.smiles)
    match = False
    for _, substructure in substructures.iterrows():
        if molecule.HasSubstructMatch(substructure.rdkit_molecule):
            matches.append(
                {
                    "chembl_id": row.molecule_chembl_id,
                    "rdkit_molecule": molecule,
                    "substructure": substructure.rdkit_molecule,
                    "substructure_name": substructure["name"],
                }
            )
            match = True
    if not match:
        clean.append(index)

matches = pd.DataFrame(matches)
egfr_data = egfr_data.loc[clean]

print(f"Numer of found unwanted structures: {len(matches)}")
print(f"Numer of compounds without unwanted sbstructures: {len(egfr_data)}")

# # Highlight structures
to_highlight = [row.rdkit_molecule.GetSubstructMatch(row.substructure) for _, row in matches.head(3).iterrows()]
Chem.Draw.MolsToGridImage(
    list(matches.head(3).rdkit_molecule),
    highlightAtomLists=to_highlight,
    legends=list(matches.head(3).substructure_name),
)

# img3 = Chem.Draw.MolsToGridImage(
#     list(matches.head(3).rdkit_molecule),
#     highlightAtomLists=to_highlight,
#     legends=list(matches.head(3).substructure_name),
# )

# img3.save("highlight_structures.png")

# # Substructure statistics: 
groups = matches.groupby("substructure_name")
group_frequencies = groups.size()
group_frequencies.sort_values(ascending=False, inplace=True)
print(group_frequencies.head(10))