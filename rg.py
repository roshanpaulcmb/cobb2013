#### IMPORTS ####
import MDAnalysis as mda
import os
from openbabel import openbabel
openbabel.obErrorLog.StopLogging()

#### FUNCTIONS ####

def calcGyrationRadius(file):
    u = mda.Universe(file)
    protein = u.select_atoms("protein") # Or specific selection
    gr = protein.radius_of_gyration() # Angstroms
    return gr

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

def convertFile(file, inputExt = "cif", outputExt = "pdb"):
    """
    Converts a CIF file to a PDB file using the core openbabel library.
    `file` must be a path to the file (not just the basename). Returns the path
    to the converted file on success or None on failure.
    """
    if not os.path.exists(file):
        print(f"Error: input file does not exist: {file}")
        return None

    obconversion = openbabel.OBConversion()
    # Set the input and output formats
    if not obconversion.SetInAndOutFormats(inputExt, outputExt):
        print("Error: Invalid input or output extension (format) specified.")
        return None

    mol = openbabel.OBMol()
    # Read the molecule from the input file
    if not obconversion.ReadFile(mol, file):
        print(f"Open Babel Error: Could not read file {file}")
        return None

    base, _ = os.path.splitext(file)
    outputFile = base + "." + outputExt
    if not obconversion.WriteFile(mol, outputFile):
        print(f"Open Babel Error: Could not write file {outputFile}")
        return None

    print(f"Successfully converted {file} to {outputFile}")
    return outputFile

#### RUN ####

def runAll(dir):
    results = {}

    for fname in sorted(os.listdir(dir)):
        name, ext = os.path.splitext(fname)
        
        if ext == ".pdb":
            path = os.path.join(dir, fname)
            try:
                gr = calcGyrationRadius(path)
                results[fname] = {"gyrationRadius": gr}
            except Exception as e:
                print(f"Failed on {fname}: {e}")
        elif ext == ".cif":
            newFname = name + ".pdb"
            path = os.path.join(dir, newFname)
            
            # If a PDB with the same basename already exists, use it and skip conversion
            if os.path.exists(path):
                print(f"Using existing PDB for {fname}: {os.path.basename(path)}")
                pass
            
            # If only the CIF exists, convert it and use it to calculate rg
            else:
                converted = convertFile(path)
                if converted is None:
                    print(f"Skipping {fname} due to conversion failure")
                    continue
                path = converted
                try:
                    gr = calcGyrationRadius(path)
                    results[os.path.basename(path)] = {"gyrationRadius": gr}
                except Exception as e:
                    print(f"Failed on {os.path.basename(path)}: {e}")    
        
        # Room to add generic else statement for converting any structure file into pdb
    
    writeCsv(os.path.join(dir, "gyrationRadii.csv"), results)
    return results

if __name__ == "__main__":
    dirs = ["ESM", "Alphafold"]
    for dir in dirs:
        runAll(dir)
    
