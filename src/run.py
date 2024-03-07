##############################################################################

# import functions
from test import *

##############################################################################
# declare path

sdf_file = "data/AcidoAcetilSalicilico.sdf"
out_file = "out.txt"

##############################################################################

# parse file
m = parse_sdf(sdf_file)

##############################################################################

# identify fragments
fragments = identify_fragments(m[0])

# print("Number of fragments:", len(fragments))
# # visualize SMILES for each fragment
# for i, fragment in enumerate(fragments, start=1):
#     print(f"Fragment {i} SMILES:", Chem.MolToSmiles(fragment))
#     print(get_rs(fragment))

##############################################################################

# get fsmiles & write to txt file
fsmiles = get_fsmiles(m[0], out_file)

# print(fsmiles)

##############################################################################

# # molecule information: atoms, bonds & connectivity

# molecule_info_list = info_sdf(sdf_file)

# # print molecule structure information 
# if molecule_info_list: 
#     first_molecule_info = molecule_info_list[0]
#     print("Atom types:", first_molecule_info['atom_types'])
#     print("Bond types:", first_molecule_info['bond_types'])
#     print("Connectivity:", first_molecule_info['connectivity'])
# else:
#     print("No molecules found in the SDF file.")

##############################################################################
