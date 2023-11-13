import sys
from rdkit import rdBase, Chem

__infile__ = "./"
__outfile__ = "./output.txt"
__mapfile__ = "./map.txt"

# Main
if __name__ == '__main__':
    args = sys.argv
    if len(args) != 4:
        print("  Invalid arguments!")
        print("    Usage: python smiles2graph.py <input_file_name> <output_file_name> <mapping_file>", file=sys.stderr)
        sys.exit()
    else:
        __infile__ = args[1]
        __outfile__ = args[2]
        __mapfile__ = args[3]

    # Read SMILES
    smilesDB = []
    f = open(__infile__, 'r')
    while True:
        line = f.readline()
        if line == "":
            break
        smilesDB.append((line.split())[1])        
    f.close()
    
    # Convert SMILES to Graph
    dump = []
    for graph_idx in range(0, len(smilesDB)):
        dump.append("t "+str(graph_idx))
        smiles = smilesDB[graph_idx]
        
        mol = Chem.MolFromSmiles(smiles)
        for atom in mol.GetAtoms():
            atom_idx = atom.GetIdx()
            atom_symbol = atom.GetSymbol()
            dump.append("v "+str(atom_idx)+" "+atom_symbol)
            
        for bond in mol.GetBonds():
            bond_idx = bond.GetIdx()
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            bond_type = bond.GetBondTypeAsDouble()
            dump.append("e "+str(bond_idx)+" "+str(begin_idx)+" "+str(end_idx)+" "+str(bond_type))

    # Output Graphs
    f = open(__outfile__, 'w')
    for idx in range(0, len(dump)):
        f.write(dump[idx]+'\n')
    f.close()

    # Output Mapping file
    f = open(__mapfile__, 'w')
    for graph_idx in range(0, len(smilesDB)):
        smiles = smilesDB[graph_idx]        
        f.write(str(graph_idx)+" "+smiles+"\n")
    f.close()
