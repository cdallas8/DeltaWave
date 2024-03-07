##############################################################################

# dependencies 
from rdkit import Chem

##############################################################################

# constants

# draft elements symbols
elements = ['C', 'N', 'O', 'S']

##############################################################################

"""
    Function to parse SDF file.

    Args:
        file_path: string - path to SDF file.

    Returns:
        list of RDKit Mol objects.
"""
def parse_sdf(file_path):
    # initialize list 
    molecules = []
    # parse file and append molecules to list 
    with open(file_path, 'rb') as f:
        m = Chem.ForwardSDMolSupplier(f)
        for mol in m:
            if mol is not None:
                molecules.append(mol)
    return molecules
    
##############################################################################

"""
    Parse an SDF file and return molecule structure information - atom & bond types, connectivity
    
    Args:
        sdf_file: string - path to SDF file
            
      Returns:
      list of dictionaries containing molecule structure information
"""
def info_sdf(file_path):
    # initialize list 
    molecules = []
    # parse file 
    with open(file_path, 'rb') as f:
        m = Chem.ForwardSDMolSupplier(f)
        for mol in m:
            if mol is not None:
                # extract atom & bond type, connectivity
                atom_types = [atom.GetSymbol() for atom in mol.GetAtoms()]
                bond_types = [bond.GetBondType() for bond in mol.GetBonds()]
                connectivity = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in mol.GetBonds()]
                # write to dict
                molecule_info = {
                    'atom_types': atom_types,
                    'bond_types': bond_types,
                    'connectivity': connectivity
                }
                molecules.append(molecule_info)
    return molecules
    
##############################################################################

"""
    Function identify fragments in a molecule

    Args:
        molecule: RDKit Mol object - molecule to cut into fragments

    Returns:
        list of RDKit Mol objects - fragments
"""
def identify_fragments(molecule):
    # get fragments 
    fragments = Chem.GetMolFrags(molecule, asMols=True)
    return fragments
 
##############################################################################

"""
    Function to obtain ring size

    Args:
        fragment: RDKit Mol object - fragment 

    Returns:
        list of int - ring size 
"""
def get_rs(fragment):
    # initialize list 
    ring_sizes = []
    
    r = Chem.GetSymmSSSR(fragment)
    for ring in r:
        ring_sizes.append(len(ring))
    return ring_sizes
    
##############################################################################

"""
    Function to generate FSMILES

    Args:
        molecule: RDKit Mol object
        fsmiles_out: string - path to output file (.txt)

    Returns:
        txt file with molecule in FSMILE format
"""
def get_fsmiles(molecule, fsmiles_out):
    # get framents
    fragments = identify_fragments(molecule)
    fsmiles = []

    # loop through all fragments
    for i, fragment in enumerate(fragments):
        # get SMILES for each fragment
        fragment_smiles = Chem.MolToSmiles(fragment)
        # get ring size
        ring_sizes = get_rs(fragment)

        # write info into FSMILES format
        fsmiles_list = []
        # check SIMLES symbols for each fragment & add to FSMILES
        for symbol in fragment_smiles:
            if symbol in elements:
                fsmiles_list.append(symbol)
                if len(ring_sizes) > 0:
                    # add ring size if present 
                    fsmiles_list.append('_' + str(ring_sizes.pop(0))) 
            elif symbol.isdigit():
                # if there is a digit, appends as last element 
                fsmiles_list[-1] += symbol
            else:
                fsmiles_list.append(symbol)
            
        #  add asterisk o denote connection points
        if i > 0:
            fsmiles_list.insert(0, '*')

        fsmiles.extend(fsmiles_list)

    fsmiles_text = ''.join(fsmiles)
    
    # write to txt file
    with open(fsmiles_out, 'w') as f:
        f.write(fsmiles_text)
        
##############################################################################