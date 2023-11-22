# GESim

## How to setup

1. Clone this repository

```bash
git clone git@github.com:LazyShion/graph_entropy.git
```

2. Install GESim

```bash
cd graph_entropy/
python3.11 -m venv .venv
source .venv/bin/activate
pip install --upgrade .
```

## How to use

```python
from rdkit import Chem
from gesim import gesim

mol1 = Chem.MolFromSmiles('CCO')
mol2 = Chem.MolFromSmiles('CCN')
mols = [mol1, mol2, Chem.MolFromSmiles('CCC')]

gesim.graph_entropy_similarity(mol1, mol2)
gesim.graph_entropy_similarity_batch(mol1, mols)
```
