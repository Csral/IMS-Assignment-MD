from rdkit import Chem
from rdkit.Chem import PDBWriter, AllChem

num_monomers = 25

smiles = "C" * (num_monomers * 2)
polymer = Chem.MolFromSmiles(smiles)
polymer = Chem.AddHs(polymer)
AllChem.EmbedMolecule(polymer)

output_file = "polyethylene.pdb"
with PDBWriter(output_file) as f:
    f.write(polymer)

print(f"Polyethylene saved to {output_file}")