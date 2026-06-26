import os
import configparser
import sys
from pathlib import Path

def replace_pdb_chain(pdb_path, new_chain):
    """
    Replaces the chain identifier of a PDB file with the provided 'new_chain',
    regardless of the original chain.
    """
    pdb_path = Path(pdb_path)
    temp_path = pdb_path.with_name(f"{pdb_path.name}.tmp")

    try:
        with open(pdb_path, 'r') as source, open(temp_path, 'w') as dest:
            for line in source:
                if line.startswith("ATOM  ") or line.startswith("HETATM"):
                    line = line[:21] + new_chain + line[22:]
                dest.write(line)

        os.replace(temp_path, pdb_path)
    except Exception as e:
        if temp_path.exists():
            temp_path.unlink(missing_ok=True)
        print(f"Error processing file {pdb_path}: {e}")

def main():
    config_file = Path('/data/config')
    base_data_dir = Path('/data')

    candidate_folders = [
        base_data_dir / 'ligands'
    ]
    folder_ligands = next((path for path in candidate_folders if path.is_dir()), None)

    config = configparser.ConfigParser()
    if not config_file.exists():
        print(f"Error: The configuration file '{config_file}' was not found.")
        sys.exit(1)
        
    config.read(config_file)
    
    try:
        pdb_receptor = config['Input_Files']['pdb_receptor']
    except KeyError:
        print("Error: The key 'pdb_receptor' is missing in config.")
        sys.exit(1)

    pdb_receptor = Path(pdb_receptor)
    if not pdb_receptor.is_absolute():
        pdb_receptor = base_data_dir / pdb_receptor

    if pdb_receptor.exists():
        if pdb_receptor.stat().st_size == 0:
            print(f"Error: The receptor file '{pdb_receptor}' is empty.")
        else:
            print(f"Changing chain of receptor '{pdb_receptor}' to 'A'...")
            replace_pdb_chain(pdb_receptor, 'A')
    else:
        print(f"Error: The receptor file '{pdb_receptor}' was not found.")


    if folder_ligands and folder_ligands.is_dir():
        pdb_files = sorted(folder_ligands.glob('*.pdb'))
        
        if not pdb_files:
            print(f"Warning: No .pdb files found in '{folder_ligands}'.")
            
        for ligand_path in pdb_files:
            if ligand_path.stat().st_size == 0:
                print(f"Error: The ligand file '{ligand_path.name}' is empty.")
                continue  
            
            print(f"Changing chain of ligand '{ligand_path.name}' to 'B'...")
            replace_pdb_chain(ligand_path, 'B')
    else:
        print("Error: No ligand folder found. Checked /data/ligands, /data/ligandos, /data/teste/ligands, /data/teste/ligandos.")
        sys.exit(1)

    print("\nPreprocessing completed successfully.")

if __name__ == "__main__":
    main()