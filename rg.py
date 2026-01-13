#### IMPORTS ####
import MDAnalysis as mda
import os
from openbabel import openbabel
openbabel.obErrorLog.StopLogging()

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

def convertFile(inputFilePath, inputFormat = "cif", outputFormat = "pdb"):
    """
    Converts a CIF file to a PDB file using the core openbabel library.
    `inputFilePath` must be a path to the file (not just the basename). Returns the path
    to the converted file on success or None on failure.
    """
    if not os.path.exists(inputFilePath):
        print(f"Error: input file does not exist: {inputFilePath}")
        return None

    obconversion = openbabel.OBConversion()
    # Set the input and output formats
    if not obconversion.SetInAndOutFormats(inputFormat, outputFormat):
        print("Error: Invalid input or output format specified.")
        return None

    mol = openbabel.OBMol()
    # Read the molecule from the input file
    if not obconversion.ReadFile(mol, inputFilePath):
        print(f"Open Babel Error: Could not read file {inputFilePath}")
        return None

    base, _ = os.path.splitext(inputFilePath)
    outputFile = base + "." + outputFormat
    if not obconversion.WriteFile(mol, outputFile):
        print(f"Open Babel Error: Could not write file {outputFile}")
        return None

    print(f"Successfully converted {inputFilePath} to {outputFile}")
    return outputFile

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
            inputPath = os.path.join(dir, fname)
            # If a PDB with the same basename already exists, use it and skip conversion
            pdbPath = os.path.splitext(inputPath)[0] + ".pdb"
            if os.path.exists(pdbPath):
                print(f"Using existing PDB for {fname}: {os.path.basename(pdbPath)}")
                pass
            
            # If only the CIF exists, convert it and use it to calculate rg
            else:
                converted = convertFile(inputPath)
                if converted is None:
                    print(f"Skipping {fname} due to conversion failure")
                    continue
                path = converted
                try:
                    rg = calcGyrationRadius(path)
                    results[os.path.basename(path)] = rg
                except Exception as e:
                    print(f"Failed on {os.path.basename(path)}: {e}")    
    
    writeGyrationRadius(os.path.join(dir, "gyrationRadii.csv"), results)
    return results

if __name__ == "__main__":
    dirs = ["ESM", "Alphafold"]
    for dir in dirs:
        runAll(dir)
    
