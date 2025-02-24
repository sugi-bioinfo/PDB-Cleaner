import os
import argparse
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, PPBuilder

def parse_arguments():
    parser = argparse.ArgumentParser(description="PDB Cleaner: Cleans and processes multiple PDB files.")
    parser.add_argument("input_folder", type=str, nargs="?", default="input_pdbs", help="Path to input PDB files (default: input_pdbs).")
    parser.add_argument("output_folder", type=str, nargs="?", default="output_pdbs", help="Path to save cleaned PDB files (default: output_pdbs).")
    parser.add_argument("-rmw", "--remove-waters", action="store_true", help="Remove waters.")
    parser.add_argument("-kh", "--keep-hydrogens", action="store_true", help="Keep hydrogens.")
    parser.add_argument("-alt", "--handle-altloc", action="store_true", help="Keep highest occupancy alternate locations.")
    parser.add_argument("-ins", "--remove-insertions", action="store_true", help="Remove insertion codes.")
    parser.add_argument("-gaps", "--report-gaps", action="store_true", help="Report sequence gaps.")
    return parser.parse_args()

def load_pdb(filepath):
    pdb_parser = PDBParser(QUIET=True)
    structure = pdb_parser.get_structure(os.path.splitext(os.path.basename(filepath))[0], filepath)
    return structure

def check_sequence_gaps(model):
    ppb = PPBuilder()
    gaps = []
    for pp in ppb.build_peptides(model, aa_only=False):
        seq = pp.get_sequence()
        for i in range(len(seq) - 1):
            if ord(seq[i + 1]) - ord(seq[i]) > 1:
                gaps.append((seq[i], seq[i + 1]))
    return gaps

def clean_pdb(structure, remove_waters, keep_hydrogens, handle_altloc, remove_insertions, report_gaps):
    model = structure[0]  # Keep first model
    cleaned_atoms = []
    sequence_gaps = []
    
    if report_gaps:
        sequence_gaps = check_sequence_gaps(model)
    
    for residue in model.get_residues():
        if residue.id[0] != " ":  # Skip HETATM records
            continue
        if remove_waters and residue.get_id()[0] == 'W':
            continue
        if residue.resname == 'MSE':  # Convert selenomethionine to methionine
            residue.resname = 'MET'
        if remove_insertions and residue.id[2] != ' ':  # Remove insertion codes
            continue
        for atom in residue:
            if atom.is_disordered():
                atom = atom.disordered_get()
                if handle_altloc and atom.get_occupancy() != max([a.get_occupancy() for a in atom.parent]):
                    continue  # Keep highest occupancy
            if atom.occupancy < 0:
                atom.set_occupancy(0.00)  # Set negative occupancy to 0.00
            else:
                atom.set_occupancy(1.00)  # Otherwise, set occupancy to 1.00
            if not keep_hydrogens and atom.element.strip() == 'H':
                continue
            cleaned_atoms.append(atom)
    return cleaned_atoms, sequence_gaps

def save_cleaned_pdb(output_folder, filename, atoms):
    os.makedirs(output_folder, exist_ok=True)
    new_filename = filename.split('.')[0][-4:] + ".pdb"  # Extract only last four characters before .pdb
    output_file = os.path.join(output_folder, new_filename)
    with open(output_file, "w") as fw:
        for i, atom in enumerate(atoms, 1):
            fw.write(f"ATOM  {i:5d} {atom.name:>4} {atom.parent.resname:>3} {atom.parent.parent.id} {atom.parent.id[1]:4d}   "
                     f"{atom.coord[0]:8.3f}{atom.coord[1]:8.3f}{atom.coord[2]:8.3f}  {atom.occupancy:6.2f} {atom.bfactor:6.2f}          {atom.element:>2}\n")

def main():
    args = parse_arguments()
    os.makedirs(args.input_folder, exist_ok=True)  # Ensure input folder exists
    os.makedirs(args.output_folder, exist_ok=True)  # Ensure output folder exists
    pdb_files = [f for f in os.listdir(args.input_folder) if f.endswith(".pdb")]
    for pdb_file in pdb_files:
        pdb_path = os.path.join(args.input_folder, pdb_file)
        structure = load_pdb(pdb_path)
        cleaned_atoms, sequence_gaps = clean_pdb(structure, args.remove_waters, args.keep_hydrogens, args.handle_altloc, args.remove_insertions, args.report_gaps)
        save_cleaned_pdb(args.output_folder, pdb_file, cleaned_atoms)
        if args.report_gaps and sequence_gaps:
            print(f"Sequence gaps in {pdb_file}: {sequence_gaps}")
        print(f"Cleaned PDB saved to {os.path.join(args.output_folder, pdb_file.split('.')[0][-4:] + '.pdb')}")

if __name__ == "__main__":
    main()
