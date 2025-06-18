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
