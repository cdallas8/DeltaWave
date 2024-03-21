##############################################################################

# Import functions
from fragmenter import *

##############################################################################
# Declare path

# Files
mol_path = "data/test_mol.sdf"
draw_path = "out/fragments.png"

# Output files 
smiles_out = "out/fragments_smiles.txt"
fsmiles_out = "out/fsmiles.txt"

##############################################################################

# Define the molecule
molecule = parse_sdf(mol_path)[0]
# molecule = Chem.MolFromSmiles('COc1ccc(-c2cccc(S(=O)(=O)N3CCN(c4nccnc4-c4ccc(OC)cc4)CC3)c2)cc1')

# Find the cuttable bonds
cuttable_bonds = find_cuttable_bonds(molecule)
print(cuttable_bonds)

# Find bonds from rings 
ring_bonds = get_ring_bond_indexes(molecule)
print(ring_bonds)

ring_bonds = merge_bonds(ring_bonds)
print(ring_bonds)

# Remove ring bonds from cuttable bonds 
bonds =  filter_bonds(cuttable_bonds, ring_bonds)
print(bonds)

# Cut the identified bonds
fragments = cut_bonds(molecule, bonds)

print("Number of fragments:", len(fragments))

for i, fragment in enumerate(fragments):
    print(f"Fragment {i+1} SMILES:", Chem.MolToSmiles(fragment))

# Write fragment SMILES to file 
fragments_file(fragments, smiles_out)

# Draw the fragments
draw_mol(fragments, draw_path)

# Generate FSMILES
fsmiles = generate_fsmiles(fragments, fsmiles_out)
print(fsmiles)

##############################################################################
