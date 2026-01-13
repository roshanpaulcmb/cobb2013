#### IMPORTS ####
import MDAnalysis as mda
import os

#### FUNCTIONS ####

def calcGyrationRadius(file):
    u = mda.Universe(file)
    protein = u.select_atoms("protein") # Or specific selection
    rg = protein.radius_of_gyration()
    print(f"Radius of Gyration: {rg} Angstroms")
    return rg

def writeGyrationRadius(fileName, results):
    with open(fileName, "w") as f:
        f.write("file, gyrationRadius\n")
        for file, gyrationRadius in results.items():
            f.write(f"{file},{gyrationRadius}\n")
    print(f"Wrote gyration radii to {fileName}")
    return None

#### RUN ####

def runAll(dir):
    results = {}

    for fname in sorted(os.listdir(dir)):
        if fname.lower().endswith(".pdb"):
            path = os.path.join(dir, fname)
            try:
                rg = calcGyrationRadius(path)
                results[fname] = rg
            except Exception as e:
                print(f"Failed on {fname}: {e}")
        elif fname.lower().endswith(".cif"):
            path = os.path.join(dir, fname)
            try:
                rg = calcGyrationRadius(path)
                results[fname] = rg
            except Exception as e:
                print(f"Failed on {fname}: {e}")    
    
    writeGyrationRadius(os.path.join(dir, "gyrationRadii.csv"), results)
    return results

if __name__ == "__main__":
    dir = "ESM"
    runAll(dir)
    
