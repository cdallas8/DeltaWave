#############################################################################

# dependencies 
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

#############################################################################

"""
    Function to parse SDF file.

    Args:
        file_path: string - path to SDF file.

    Returns:
        list of RDKit Mol objects - list of molecules parsed from the SDF file.
"""
def parse_sdf(file_path):
    # Initialize list 
    molecules = []
    # Parse file and append molecules to list 
    with open(file_path, 'rb') as f:
        m = Chem.ForwardSDMolSupplier(f)
        for mol in m:
            if mol is not None:
                molecules.append(mol)
    return molecules
    

"""
    Function to find cuttable bonds in a molecule.

    Args:
        molecule: RDKit Mol object. 

    Returns:
        indices of cuttable bonds.
"""
def find_cuttable_bonds(molecule):
    # Get the ring information for the molecule
    ri = molecule.GetRingInfo()
    
    # Get number of rings 
    num_rings = ri.NumRings()
    print("Number of rings in the molecule:", num_rings)

    # Initialize a list to store cuttable bonds
    cuttable_bonds = []

    # Iterate through the bonds in the molecule
    for bond in molecule.GetBonds():
        # Check if the bond is a single bond
        if bond.GetBondType() == Chem.BondType.SINGLE:
            # Get the atoms at both ends of the bond
            begin_atom = molecule.GetAtomWithIdx(bond.GetBeginAtomIdx())
            end_atom = molecule.GetAtomWithIdx(bond.GetEndAtomIdx())
            
            # Get index of the bond
            bondx = bond.GetIdx()
            print("idx", bondx)

            # Check if both ends of the bond are not connected to hydrogen atoms
            if not (begin_atom.GetAtomicNum() == 1 or end_atom.GetAtomicNum() == 1):
                # Check if at least one end of the bond is attached to a ring
                if (ri.NumAtomRings(begin_atom.GetIdx()) > 0 or ri.NumAtomRings(end_atom.GetIdx()) > 0):
                    # # Check if at least one end of the bond is not part of the ring 
                    # if (ri.NumAtomRings(begin_atom.GetIdx()) == 0 or ri.NumAtomRings(end_atom.GetIdx()) == 0) or (begin_atom.IsInRing() or end_atom.IsInRing()):
                    
                        # Append the bond index to the list of cuttable bonds
                        cuttable_bonds.append(bond.GetIdx())

    return cuttable_bonds
    
#############################################################################

"""
    Function to fragment molecule given bonds indices.

    Args:
        molecule: RDKit Mol object. 
        bond_indices: list of bond indices specifying which bonds to cut for fragmentation.

    Returns:
        list of RDKit Mol objects - the fragments resulting from cutting the specified bonds.
"""
def cut_bonds(molecule, bond_indices):
    # Fragment the molecule on the specified bonds
    fragments = Chem.FragmentOnBonds(molecule, bond_indices)
    frags = Chem.GetMolFrags(fragments, asMols = True)
    return frags

#############################################################################

"""
    Function to draw fragments of a molecule.

    Args:
        fragments: list of RDKit Mol objects - fragments of the molecule to be drawn.
        out_path: str - path to save the .png file.
    
    Returns:
        None
"""
def draw_mol(fragments, out_path):
  Draw.MolsToGridImage(fragments).save(out_path)
  
#############################################################################

"""
Function to get the bond indexes for each ring in a molecule.

Args:
    molecule: RDKit Mol object. 

Returns:
    list of lists - each inner list contains the bond indexes for a ring in the molecule.
"""
def get_ring_bond_indexes(molecule):
    # Get the ring information for the molecule
    ri = molecule.GetRingInfo()
    
    # Create a list to store all bond indexes in each ring
    all_ring_bond_indexes = []

    # Iterate through the rings in the molecule
    for ring in ri.AtomRings():
        # Create a list to store bond indexes in the current ring
        bond_indexes = []
        # Iterate through the atoms in the ring
        for i in range(len(ring)):
            # Get the bond between the current atom and the next atom in the ring
            bond = molecule.GetBondBetweenAtoms(ring[i], ring[(i + 1) % len(ring)])
            # Append the bond index to the list
            bond_indexes.append(bond.GetIdx())
        # Append the list of bond indexes for the current ring to the list of all ring bond indexes
        all_ring_bond_indexes.append(bond_indexes)

    return all_ring_bond_indexes
    
#############################################################################

"""
Function to merge bond indexes from all rings into a single list.

Args:
    all_ring_bond_indexes: list of lists. 

Returns:
    list - merged list containing all bond indexes from all rings.
"""
def merge_bonds(all_ring_bond_indexes):
    merged_idx = []
    for bond_indexes in all_ring_bond_indexes:
        merged_idx.extend(bond_indexes)
    
    return merged_idx
    
#############################################################################

"""
Function to filter cuttable bonds based on ring bond indices.

Args:
    cuttable_bonds: list - list of bond indices that are potentially cuttable.
    ring_bond_indexes: list - list of bond indices belonging to rings.

Returns:
    list - filtered list of cuttable bond indices that are not part of any ring.
"""
def filter_bonds(cuttable_bonds, ring_bond_indexes):
    # Initialize a list to store the filtered cuttable bonds
    filtered_cuttable_bonds = []
    
    # Iterate through the cuttable bonds
    for bond_index in cuttable_bonds:
        # Check if the bond index is not in the list of ring bond indexes
        if bond_index not in ring_bond_indexes:
            # Add the bond index to the filtered list
            filtered_cuttable_bonds.append(bond_index)
    
    return filtered_cuttable_bonds

#############################################################################

"""
Function to write fragments to a file in SMILES format.

Args:
    fragments: list - list of RDKit Mol objects representing fragments.
    filename: str - path to the output file.

Returns:
    None
""""
def fragments_file(fragments, filename):
    with open(filename, 'w') as f:
        for i, fragment in enumerate(fragments):
            # Convert the fragment to SMILES format
            fragment_smiles = Chem.MolToSmiles(fragment)
            # Write the fragment SMILES to the file
            f.write(fragment_smiles)
            # Write the separator 'sep' after each fragment except for the last one
            if i < len(fragments) :
                f.write("'sep'")

#############################################################################

"""
Function to get the size of rings in a fragment.

Args:
    fragment: RDKit Mol object - fragment for which ring sizes are to be obtained.

Returns:
    int - size of the largest ring in the fragment, or 0 if the fragment contains no rings.
"""
def get_ringsize(fragment):
    # Get the ring information for the molecule
    ri = fragment.GetRingInfo()
    # Check if the fragment contains rings
    if ri.NumRings() > 0:
        # Extract ring sizes from the ring information
        for ring in ri.AtomRings():
            ring_size = len(ring)
    else: 
        ring_size = 0
    return ring_size
    
#############################################################################

"""
Function to generate FSMILES for fragments.

Args:
    fragments: list - list of RDKit Mol objects representing fragments.

Returns:
    str - FSMILES for the fragments.
"""
def generate_fsmiles(fragments):
    fsmiles_list= []
    for fragment in fragments:
        s = get_ringsize(fragment)
        f = Chem.MolToSmiles(fragment)
  
        updated_fragment = ''
        inside_bracket = False
        for char in f:
            if char == '[':
                inside_bracket = True
                updated_fragment += char
            elif char == ']':
                inside_bracket = False
                updated_fragment += char + '_0'
            elif not inside_bracket and char.isalpha():
                updated_fragment += char + '_' + str(s)
            else:
                updated_fragment += char

        fsmiles_list.append(updated_fragment)
    
    # format 
    fsmiles_list = "'sep_0'".join(fsmiles_list)     
    
    return fsmiles_list
    
#############################################################################
