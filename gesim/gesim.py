import os
from typing import List
import tempfile

from rdkit import Chem
from gesim import convert, gdb, graph_entropy


def graph_entropy_similarity(
    mol1: Chem.Mol,
    mol2: Chem.Mol,
    r: int=4,
    return_graph_entropy: bool=False) -> float:
    """Calculate the graph entropy similarity between two molecules.

    Args:
        mol1 (Chem.Mol): Mol object as defined by RDKit
        mol2 (Chem.Mol): Mol object as defined by RDKit
        r (int, optional): Radius parameter of fingerprint calculation.
                           Defaults to 4.
        return_graph_entropy (bool, optional): If set to `True`, the function returns
                                               the raw graph entropy. If `False`, it 
                                               returns the similarity measure, which 
                                               is 1 minus the graph entropy.
                                               Defaults to False.

    Returns:
        float: the graph entropy similarity measure if `return_graph_entropy` is `False`.
               Otherwise, returns the graph entropy value itself.
    """

    mol_list = [mol1, mol2]
    with tempfile.TemporaryDirectory(prefix='gesim_') as tmpdir:
        tmp_file = os.path.join(tmpdir, 'graph.txt')
        tmp_bin_file = os.path.join(tmpdir, 'graph.bin')
        create_graph_from_mols(mol_list, tmp_file)
        convert.convert_graph_to_binary(tmp_file, tmp_bin_file)
        db = gdb.GraphDB(tmp_bin_file)
    ge_calculator = graph_entropy.GraphEntropy(db, r)
    ge_value = ge_calculator.graph_entropy(0, 1)
    return ge_value if return_graph_entropy else 1 - ge_value
        

def graph_entropy_similarity_batch(
    mol_query: Chem.Mol,
    mols: List[Chem.Mol],
    r: int=4,
    return_graph_entropy: bool=False) -> List[float]:
    """Calculate the graph entropy similarity for a query molecule against a batch of
       molecules.

    Args:
        mol_query (Chem.Mol): Query RDKit Mol object to be compared against the list
        mols (List[Chem.Mol]): List of RDKit Mol objects to compare with the query molcule.
        r (int, optional): Radius parameter of fingerprint calculation.
                           Defaults to 4.
        return_graph_entropy (bool, optional): If set to `True`, the function returns 
                                               the list of raw graph entropy values. 
                                               If `False`, it returns the list of 
                                               the similarity measures, which is 1 
                                               minus the graph entropy values.
                                               Defaults to False.

    Returns:
        List[float]: List of similarity measures between the query molecule and each
                     molecule in the list if `return_graph_entropy` is `False`. If 
                     `True`, returns a list of graph entropy values for each comparison.
    """
    mol_list = [mol_query] + mols
    with tempfile.TemporaryDirectory(prefix='gesim_') as tmpdir:
        tmp_file = os.path.join(tmpdir, 'graph.txt')
        tmp_bin_file = os.path.join(tmpdir, 'graph.bin')
        create_graph_from_mols(mol_list, tmp_file)
        convert.convert_graph_to_binary(tmp_file, tmp_bin_file)
        db = gdb.GraphDB(tmp_bin_file)
    ge_calculator = graph_entropy.GraphEntropy(db, r)
    ge_values = ge_calculator.graph_entropy_all(0)
    return ge_values if return_graph_entropy else [1 - x for x in ge_values]


def create_graph_from_mols(
    mol_list: List[Chem.Mol],
    output_graph_file: str) -> None:
    """Generate graph representations of a given list of molecules and save them to a file.

    Args:
        mol_list (List[Chem.Mol]): List of RDKit Mol objects to be converted into graph
                                   representations.
        output_graph_file (str): Path to output of graph data
    """
    output_list = []
    for i, mol in enumerate(mol_list):
        crs_graph = ""
        crs_graph += f"t {i}\n"
        for a in mol.GetAtoms():
            crs_graph += f"v {a.GetIdx()} {a.GetSymbol()}\n"
        for b in mol.GetBonds():
            crs_graph += f"e {b.GetIdx()} {b.GetBeginAtomIdx()} {b.GetEndAtomIdx()} {b.GetBondTypeAsDouble()}\n"
        output_list.append(crs_graph)

    with open(output_graph_file, 'w') as f:
        f.writelines(output_list)
 