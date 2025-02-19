from rdkit import Chem
from rdkit.Chem import AllChem

suppl = Chem.SDMolSupplier('teflon.sdf')

for i, mol in enumerate(suppl):
    if mol is not None:
        if mol.GetNumConformers() == 0:
            print(f"Molecule {i+1} has no 3D coordinates. Generating 3D coordinates...")
            mol = Chem.AddHs(mol)  
            AllChem.EmbedMolecule(mol)
            AllChem.UFFOptimizeMolecule(mol)

        output_pdb = f'molecule_{i+1}.pdb'
        Chem.MolToPDBFile(mol, output_pdb)
        print(f"Saved molecule {i+1} to {output_pdb}")
    else:
        print(f"Failed to process molecule {i+1}")