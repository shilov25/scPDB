from Bio import PDB
import os
import sys
import numpy as np

def split_complex_into_binders(pdb_file, chains, output_dir, run_name):
    """
    Splits PDB input file into two binders based on provided chain map. 
    Example: chains = "AB_CD" means chains A and B are binder 1, chains C and D are binder 2. 
    Saves two new PDB files in the output directory with results (this can be turned off).
    Returns coordinates for binder 1 and binder 2 as numpy arrays.
    """
    #parse chain map
    binder1_chains, binder2_chains = chains.split("_")

    #create new PDB files for each binder, these are necessary for PyMOL mesh generation
    binder1_pdb = os.path.join(output_dir, run_name, f"{run_name}_{binder1_chains}.pdb")
    binder2_pdb = os.path.join(output_dir, run_name, f"{run_name}_{binder2_chains}.pdb")

    binder1_coords = []
    binder2_coords = []
    binder1_lines = []
    binder2_lines = []
    
    with open(pdb_file, 'r') as f:
        for line in f:
            #convert MSE to MET
            if line[:6] == "HETATM" and line[17:20] == "MSE":
                line = line.replace("HETATM", "ATOM  ")
                line = line.replace("MSE", "MET")

            #collect chain information and coordinates
            if line[:4] == "ATOM":
                chain = line[21:22]
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                coord = np.array([x, y, z])

                #match chain to binder 1 or binder 2
                if chain in binder1_chains:
                    binder1_coords.append(coord)
                elif chain in binder2_chains:
                    binder2_coords.append(coord)
                
                if chain in binder1_chains:
                    binder1_lines.append(line)
                elif chain in binder2_chains:
                    binder2_lines.append(line)

            elif line[:3] == "TER":
                if chain in binder1_chains:
                    binder1_lines.append(line)
                if chain in binder2_chains:    
                    binder2_lines.append(line)

    binder1_coords = np.array(binder1_coords)
    binder2_coords = np.array(binder2_coords)

    with open(binder1_pdb, 'w') as f:
        for line in binder1_lines:
            f.write(line)
        f.write("END\n")
            
    with open(binder2_pdb, 'w') as f:
        for line in binder2_lines:
            f.write(line)
        f.write("END\n")

    return binder1_chains, binder2_chains, binder1_coords, binder2_coords, binder1_pdb, binder2_pdb
