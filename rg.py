import MDAnalysis as mda
u = mda.Universe("ESM/random1.pdb")
protein = u.select_atoms("protein") # Or specific selection
rg = protein.radius_of_gyration()
print(f"Radius of Gyration: {rg} Angstroms")

