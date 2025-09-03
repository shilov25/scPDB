import argparse
import os
import pandas as pd
from scPDB.pdb_parser import split_complex_into_binders
from scPDB.shape_complementarity import calculate_sc

def run(pdb_file, chains, output_dir, run_name, w_constant, surface_point_density, interface_thresh, save_obj_files, save_binder_files):
    binder1_chains, binder2_chains, _, _, binder1_pdb, binder2_pdb = split_complex_into_binders(
        pdb_file, 
        chains, 
        output_dir, 
        run_name,
        )
    binder1_to_binder2_sc, binder2_to_binder1_sc, sc, b1_interface_area, b2_interface_area = calculate_sc(
        binder1_chains,
        binder2_chains,
        binder1_pdb,
        binder2_pdb,
        run_name,
        output_dir, 
        interface_thresh, 
        surface_point_density, 
        w_constant,
        save_obj_files,
        save_binder_files,
        )
    
    binder1_chains, binder2_chains = chains.split("_")

    shape_complementarity_data = {
        "PDB_ID": [os.path.splitext(os.path.basename(pdb_file))[0]],
        "Binder 1 Chains": [binder1_chains],
        "Binder 2 Chains": [binder2_chains],
        "Binder 1 Sc": [round(binder1_to_binder2_sc, 3)],
        "Binder 1 Interface Area": [round(b1_interface_area)],
        "Binder 2 Sc": [round(binder2_to_binder1_sc, 3)],
        "Binder 2 Interface Area": [round(b2_interface_area)],
        "Complex Sc": [round(sc, 3)]
    }

    sc_df = pd.DataFrame(shape_complementarity_data)
    return sc_df

def batch_run(batch_csv, output_dir, w_constant, surface_point_density, interface_thresh, save_obj_files, save_binder_files):
    batch_df = pd.read_csv(batch_csv)
    #initialize run_name column
    batch_df['run_name'] = batch_df['pdb_file'].apply(lambda x: os.path.splitext(os.path.basename(x))[0])

    all_scores = []
    for i, row in batch_df.iterrows():
        pdb_file = row['pdb_file']
        chains = row['chains']
        run_name = row['run_name']

        print(f"Processing {run_name}.")

        if not os.path.exists(pdb_file):
            print(f"{run_name} PDB file not found.")
            continue
        if "_" not in chains:
            raise ValueError(f"{run_name} has invalid chain format. Valid example: AB_CD")
        os.makedirs(os.path.join(output_dir, run_name), exist_ok=True)

        sc_df = run(pdb_file, chains, output_dir, run_name, w_constant, surface_point_density, interface_thresh, save_obj_files, save_binder_files)
        output_file = os.path.join(output_dir, run_name, f"{run_name}_sc_scores.csv")
        sc_df.to_csv(output_file, index=False)

        all_scores.append(sc_df)
    
    combined_sc_df = pd.concat(all_scores, ignore_index=True)
    batch_csv_basename = os.path.splitext(os.path.basename(batch_csv))[0]
    combined_output_file = os.path.join(output_dir, f"{batch_csv_basename}_sc_scores.csv")
    print(f"Saved batch sc score results to {combined_output_file}")
    combined_sc_df.to_csv(combined_output_file, index=False)

    return combined_sc_df

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    #basic arguments
    argparser.add_argument("--run_name", type=str, default=None, help="Naming convention for output files. If not provided, uses PDB filename.")
    argparser.add_argument("--output_dir", type=str, default="./", help="Path to output directory.")
    argparser.add_argument("--batch_csv", type=str, help="Path to .csv file containing 'pdb_file' and 'chains' columns for batch processing.")

    #arguments relevant to shape complementarity calculation
    argparser.add_argument("--pdb_file", type=str, help="Path to PDB file.")
    argparser.add_argument("--chains", type=str, help="Chains for which to calculate shape complementary score.")
    argparser.add_argument("--w_constant", type=float, default=0.5, help="Non-zero scalar constant for defining rate of sc score decay.")
    argparser.add_argument("--surface_point_density", type=int, default=30, help="Number of points/Angstrom^2 to sample.")
    argparser.add_argument("--interface_thresh", type=float, default=1.5, help="Distance threshold (in Angstroms) for defining interaction interface between two binders.")
    argparser.add_argument("--save_obj_files", action="store_true", help="If provided, saves PyMOL interface mesh outputs as individual .obj files.")
    argparser.add_argument("--save_binder_files", action="store_true", help="If provided, saves Binder 1 and Binder 2 as individual .pdb files.")
    args = argparser.parse_args()

    #for single PDB processing
    if args.pdb_file:
        if os.path.exists(args.pdb_file):
            PDB_ID = os.path.splitext(os.path.basename(args.pdb_file))[0]
        if not args.run_name:
            args.run_name = PDB_ID
        else:
            raise ValueError("PDB file does not exist.")
        
        if not args.chains:
            raise ValueError("Chains not provided.")
        if args.chains and "_" not in args.chains:
            raise ValueError("Chains must be provided in the format 'AB_CD' where AB are chains for Binder 1 and CD are chains for Binder 2.")
        
        if args.interface_thresh <= 0:
            raise ValueError("Interface threshold must be > 0 Angstroms.")
        if args.surface_point_density <= 0:
            raise ValueError("At least 1/Angstrom^2 surface points are required to calculate shape complementarity.")
        if args.w_constant <= 0:
            raise ValueError("w must be > 0")
        
        if os.path.exists(args.output_dir):
            os.makedirs(os.path.join(args.output_dir, args.run_name), exist_ok=True)
        else:
            raise ValueError("Output directory path does not exist.")
        
        sc_df = run(
        args.pdb_file, 
        args.chains, 
        args.output_dir, 
        args.run_name, 
        args.w_constant, 
        args.surface_point_density, 
        args.interface_thresh, 
        args.save_obj_files,
        args.save_binder_files,
        )

        sc_score_file = os.path.join(args.output_dir, args.run_name, f"{args.run_name}_sc_scores.csv")
        sc_df.to_csv(sc_score_file, index=False)
        print(f"Saved shape complementarity scores for {PDB_ID} to {sc_score_file}.")

    #for batch PDB processing
    elif args.batch_csv:
        if not os.path.exists(args.batch_csv):
            raise ValueError("CSV file does not exist.")

        if args.interface_thresh <= 0:
            raise ValueError("Interface threshold must be > 0 Angstroms.")
        if args.surface_point_density <= 0:
            raise ValueError("At least 1/Angstrom^2 surface points are required to calculate shape complementarity.")
        if args.w_constant <= 0:
            raise ValueError("w must be > 0")
    
        if not os.path.exists(args.output_dir):
            raise ValueError("Output directory path does not exist.")

        combined_sc_df = batch_run(
            args.batch_csv, 
            args.output_dir, 
            args.w_constant, 
            args.surface_point_density, 
            args.interface_thresh, 
            args.save_obj_files,
            args.save_binder_files,
        )
            