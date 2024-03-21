#############################################################################

# dependencies 
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

#############################################################################

def find_cuttable_bonds(molecule):
    # Get the ring information for the molecule
    ri = molecule.GetRingInfo()
    
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
                    # Check if at least one end of the bond is not part of the ring 
                    if (ri.NumAtomRings(begin_atom.GetIdx()) == 0 or ri.NumAtomRings(end_atom.GetIdx()) == 0) or (begin_atom.IsInRing() or end_atom.IsInRing()):
                    
                    # Append the bond index to the list of cuttable bonds
                        cuttable_bonds.append(bond.GetIdx())

    return cuttable_bonds
    
#############################################################################

def cut_bonds(molecule, bond_indices):
    # Fragment the molecule on the specified bonds
    fragments = Chem.FragmentOnBonds(molecule, bond_indices)
    frags = Chem.GetMolFrags(fragments, asMols = True)
    return frags

#############################################################################

def draw_mol(fragments):
  Draw.MolsToGridImage(fragments).save("fragments.png")
  
#############################################################################

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

def merge_bonds(all_ring_bond_indexes):
    merged_idx = []
    for bond_indexes in all_ring_bond_indexes:
        merged_idx.extend(bond_indexes)
    
    return merged_idx
    
#############################################################################

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
