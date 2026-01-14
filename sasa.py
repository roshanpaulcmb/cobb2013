#### IMPORTS ####
import freesasa
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
import os

#### FUNCTIONS ####

def calcSasa(file, method = "LeeRichards"):
    if method == "LeeRichards":
        structure = freesasa.Structure(file)
        totalSasa = freesasa.calc(structure)
        sasa = freesasa.classifyResults(totalSasa, structure)
    
    elif method == "ShrakeRupley":
        params = freesasa.Parameters({
            'algorithm': freesasa.ShrakeRupley,
            'n-points': 500  # default resolution is 100 for both freesasa and biopython, increasing this improves accuracy
        })
        params.setProbeRadius(1.4) # default probe radius is 1.4 A because we are rolling water
        structure = freesasa.Structure(file)
        totalSasa = freesasa.calc(structure, params)
        sasa = freesasa.classifyResults(totalSasa, structure)    
    
    elif method == "Bio.PDB":
        parser = PDBParser()
        base, _ = os.path.splitext(file)
        structure = parser.get_structure(base, file)
        total, sasa = freesasa.calcBioPDB(structure)
    
    else:
        print("Method must be 'LeeRichards', 'ShrakeRupley', or 'Bio.PDB'")

    return sasa

#### SCRIPTS ####
# Biopython ShrakeRupley

# file = "ESM/random1.pdb" 

# p = PDBParser()
# pdb = p.get_structure("random0","ESM/random1.pdb")
# sr = ShrakeRupley()
# sr.compute(pdb, level="S")
# biopy_total = pdb.sasa
# print("biopython total sasa using ShrakeRupley",biopy_total)

# # freesasa Lee-Richards
# print ("freesasa results:")
# structure = freesasa.Structure("ESM/random1.pdb")
# result = freesasa.calc(structure)
# area_classes = freesasa.classifyResults(result, structure)

# print("Total : ", result.totalArea(), "using Lee-Richards")
# for key in area_classes:
#     print(key, ": ",  area_classes[key], "using Lee-Richards")

# # freesasa Shrake-Rupley
# params_sr = freesasa.Parameters({
#     'algorithm': freesasa.ShrakeRupley,
#     'n-points': 500  # default resolution is 100 for both freesasa and biopython, increasing this improves accuracy
# })
# params_sr.setProbeRadius(1.4)  # default probe radius is 1.4 A because we assume we are rolling water
# structure_sr = freesasa.Structure("ESM/random1.pdb")
# result_sr = freesasa.calc(structure_sr, params_sr)
# area_classes_sr = freesasa.classifyResults(result_sr, structure_sr)    
# print("Total : ", result_sr.totalArea(), "using ShrakeRupley")    
# for key in area_classes_sr:
#     print(key, ": ", area_classes_sr[key], "using ShrakeRupley")      

def writeCsv(fileName, results):
    # Collect headers from the first entry
    features = next(iter(results.values()))
    headers = ["file"] + list(features.keys())

    with open(fileName, "w") as f:
        # Write header
        f.write(", ".join(headers) + "\n")

        # Write rows
        for file, features in results.items():
            values = [str(features[h]) for h in features.keys()]
            f.write(f"{file}, " + ", ".join(values) + "\n")

    print(f"Wrote {fileName}")
    return None

#### RUN ####

def runAll(dir):
    results = {}
    
    for fname in sorted(os.listdir(dir)):
        name, ext = os.path.splitext(fname)
        
        if ext == ".pdb":
            
            path = os.path.join(dir, fname)
            sasa = calcSasa(path, method = "LeeRichards")
            results[fname] = sasa
            # calcSasa(path, method = "ShrakeRupley")
            # calcSasa(path, method = "Bio.PDB")
    writeCsv(os.path.join(dir, "sasaLeeRichards.csv"), results)
    return None

if __name__ == "__main__":
    dirs = ["ESM", "Alphafold"]
    for dir in dirs:
        runAll(dir)