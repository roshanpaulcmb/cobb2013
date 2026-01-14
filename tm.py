import numpy as np
from Bio.PDB import PDBParser
from tmtools.io import get_structure, get_residue_data
from tmtools import tm_align
from tmtools.testing import get_pdb_path
p = PDBParser()
pdb_esm = p.get_structure("random9","ESM/random9.pdb")
chain_esm = next(pdb_esm.get_chains())

pdb_af3 = p.get_structure("fold_random9_model_0","AlphaFold/fold_random9_model_0.pdb")
chain_af3 = next(pdb_af3.get_chains())

coords_esm, seq_esm = get_residue_data(chain_esm)
coords_af3, seq_af3 = get_residue_data(chain_af3)
print("ESM coords shape:", coords_esm.shape)
print("AF3 coords shape:", coords_af3.shape)
res = tm_align(coords_esm, coords_af3, seq_esm, seq_af3)
print(res.tm_norm_chain1)
print("RMSD:", res.rmsd)
