##############################################################################

# import functions
from fragmenter import *

##############################################################################
# declare path


##############################################################################

# Define the molecule
molecule = Chem.MolFromSmiles('COc1ccc(-c2cccc(S(=O)(=O)N3CCN(c4nccnc4-c4ccc(OC)cc4)CC3)c2)cc1')

# Find the cuttable bonds
cuttable_bonds = find_cuttable_bonds(molecule)
print(cuttable_bonds)

ring_bonds = get_ring_bond_indexes(molecule)
print(ring_bonds)

ring_bonds = merge_bonds(ring_bonds)
print(ring_bonds)

bonds =  filter_bonds(cuttable_bonds, ring_bonds)
print(bonds)

# Cut the identified bonds
fragments = cut_bonds(molecule, bonds)

print("Number of fragments:", len(fragments))

for i, fragment in enumerate(fragments):
    print(f"Fragment {i+1} SMILES:", Chem.MolToSmiles(fragment))

draw_mol(fragments)

##############################################################################
