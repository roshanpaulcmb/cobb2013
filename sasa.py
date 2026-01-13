import freesasa
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

# Biopython ShrakeRupley
p = PDBParser()
pdb = p.get_structure("random0","ESM/random1.pdb")
sr = ShrakeRupley()
sr.compute(pdb, level="S")
biopy_total = pdb.sasa
print("biopython total sasa using ShrakeRupley",biopy_total)

# freesasa Lee-Richards
print ("freesasa results:")
structure = freesasa.Structure("ESM/random1.pdb")
result = freesasa.calc(structure)
area_classes = freesasa.classifyResults(result, structure)

print("Total : ", result.totalArea(), "using Lee-Richards")
for key in area_classes:
    print(key, ": ",  area_classes[key], "using Lee-Richards")

# freesasa Shrake-Rupley
params_sr = freesasa.Parameters({
    'algorithm': freesasa.ShrakeRupley,
    'n-points': 500  # default resolution is 100 for both freesasa and biopython, increasing this improves accuracy
})
params_sr.setProbeRadius(1.4)  # default probe radius is 1.4 A because we assume we are rolling water
structure_sr = freesasa.Structure("ESM/random1.pdb")
result_sr = freesasa.calc(structure_sr, params_sr)
area_classes_sr = freesasa.classifyResults(result_sr, structure_sr)    
print("Total : ", result_sr.totalArea(), "using ShrakeRupley")    
for key in area_classes_sr:
    print(key, ": ", area_classes_sr[key], "using ShrakeRupley")      
