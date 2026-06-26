import configparser
import subprocess
import sys
import csv
import shutil
import time
import os
import random
import hashlib
import pandas as pd # type: ignore

from pathlib import Path

import script_mat


# ============================================================
# CONSTANTS
# ============================================================

BASE_DATA = Path("/data")
HOST_DATA = Path(os.environ.get("HOST_DATA_DIR", "/data"))

SCRIPT_PRE = Path(__file__).parent / "script_pre.py"
LOG_FILE = BASE_DATA / "results.csv"

ROSETTA_IMAGE = "pegi3s/rosettadock"
CCP4_IMAGE = "pegi3s/ccp4"
MAX_PISA_PDB_NAME_LEN = len("production_6194_AF_t2_0002.pdb")


# ============================================================
# UTILITIES
# ============================================================

def run_command(cmd, description):
    """Run subprocess and log errors, but DO NOT sys.exit(1)."""
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"\n[ERROR] {description} failed")
        if result.stderr:
            print(result.stderr.strip())

    return result


def rename_rosetta_output(prefix, new_name, output_dir=BASE_DATA, search_dir=BASE_DATA):
    """Rename Rosetta output to predictable filename. Returns None if fail."""
    candidates = list(search_dir.glob(f"{prefix}*.pdb"))

    if not candidates:
        print(f"Rosetta produced no PDB")
        return None

    latest = max(candidates, key=lambda p: p.stat().st_mtime)

    output_dir.mkdir(parents=True, exist_ok=True)
    new_path = output_dir / new_name
    latest.rename(new_path)

    return new_path


def ensure_log(region_ids):
    if not LOG_FILE.exists():
        with open(LOG_FILE, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(
                ["Ligand","Iteration","Status","Iteration_Time(s)",
                 "Metric","TSum","FSum"]
                + [f"Perc{r}" for r in region_ids]
            )


def log_result(ligand, iteration, status, elapsed, data=None, region_ids=None):
    region_ids = region_ids or []
    with open(LOG_FILE, "a", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        if data is not None:
            writer.writerow(
                [ligand, iteration, status, elapsed,
                 round(data.get("Metrica", 0), 4),
                 round(data.get("TSum", 0), 4),
                 round(data.get("FSum", 0), 4)]
                + [round(data.get(f"P{r}", 0), 2) for r in region_ids]
            )
        else:
            writer.writerow([ligand, iteration, status, elapsed] + ["N/A"] * (3 + len(region_ids)))


def build_pisa_input_name(pdb_name):
    """Return a CCP4/PISA-safe input filename capped to a known-safe length."""
    if len(pdb_name) <= MAX_PISA_PDB_NAME_LEN:
        return pdb_name

    original = Path(pdb_name)
    digest = hashlib.sha1(pdb_name.encode("utf-8")).hexdigest()[:8]
    stem = original.stem

    # Preserve Rosetta model index when present (e.g. _0001).
    model_suffix = ""
    last_token = stem.rsplit("_", 1)[-1]
    if len(last_token) == 4 and last_token.isdigit():
        model_suffix = f"_{last_token}"

    if stem.startswith("production_"):
        safe_name = f"production_{digest}{model_suffix}.pdb"
    else:
        safe_name = f"pisa_{digest}{model_suffix}.pdb"

    if len(safe_name) > MAX_PISA_PDB_NAME_LEN:
        # Final hard cap fallback, preserving .pdb extension.
        safe_name = f"{safe_name[:MAX_PISA_PDB_NAME_LEN - 4]}.pdb"

    return safe_name


def evaluate_pisa_outputs(folder_out, config_path):
    """Evaluate PISA TSV outputs and return first valid data payload."""
    tsv_candidates = sorted(folder_out.rglob("interface_*_info.tsv"))

    for tsv_file in tsv_candidates:
        try:
            success, data = script_mat.evaluate(tsv_file, config_path)
        except Exception as e:
            print(f"[ERROR] Math evaluation error for {tsv_file}: {e}")
            continue

        if data is not None:
            return success, data, tsv_file, tsv_candidates

    return False, None, None, tsv_candidates

def extract_structure_summaries(tsv_file, output_dir):
    """Extract Structure1 (lines 8-11) and Structure2 (lines 13-16) from interface TSV."""
    try:
        with open(tsv_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()

        # Structure1: lines 8 to 11 (index 7 to 10)
        if len(lines) >= 11:
            with open(output_dir / "Structure1_summary.tsv", 'w', encoding='utf-8') as f:
                f.writelines(lines[7:11])

        # Structure2: lines 13 to 16 (index 12 to 15)
        if len(lines) >= 16:
            with open(output_dir / "Structure2_summary.tsv", 'w', encoding='utf-8') as f:
                f.writelines(lines[12:16])

    except Exception as e:
        print(f"[ERROR] Could not extract structure summaries: {e}")

def ensure_control_csv(base_dir):
    """If control file is .xlsx, convert to .csv. Returns csv path or None."""
    csv_file = base_dir / "control_results.csv"
    xlsx_file = base_dir / "control_results.xlsx"

    if csv_file.exists():
        return csv_file

    if xlsx_file.exists():
        try:
            df = pd.read_excel(xlsx_file)
            df.to_csv(csv_file, index=False)
            print(f"Converted {xlsx_file.name} to CSV.")
            return csv_file
        except Exception as e:
            print(f"[WARNING] Could not convert {xlsx_file.name}: {e}")
            return None

    return None

def generate_summary(control_file, results_file, output_file, region_ids):
    """Read GeneIDs from control file, find best result for each in results.csv."""

    try:
        control_df = pd.read_csv(control_file)
        gene_ids = control_df.iloc[:, 0].astype(str).tolist()
    except Exception as e:
        print(f"[WARNING] Could not read control file '{control_file}': {e}")
        return

    if not results_file.exists():
        print(f"[WARNING] Results file '{results_file}' not found. Skipping summary.")
        return

    try:
        results_df = pd.read_csv(results_file)
    except Exception as e:
        print(f"[WARNING] Could not read results file: {e}")
        return

    def normalize(s):
        return str(s).replace('_', '').lower()

    results_df['_norm'] = results_df['Ligand'].apply(lambda x: normalize(x))
    sorted_regions = sorted(region_ids)

    rows = []
    for gene_id in gene_ids:
        norm_id = normalize(gene_id)
        matches = results_df[results_df['_norm'].str.contains(norm_id, na=False)]

        if matches.empty:
            continue

        # Prefer Approved with lowest Metric
        approved = matches[matches['Status'] == 'Approved']
        if not approved.empty:
            best = approved.loc[approved['Metric'].idxmin()]
        else:
            valid = matches[matches['Metric'].notna() & (matches['Metric'] > 0)]
            if valid.empty:
                continue
            best = valid.loc[valid['Metric'].idxmin()]

        row = {'GeneID': gene_id}
        for i, rid in enumerate(sorted_regions, 1):
            col = f'P{rid}'
            row[f'region{i}'] = round(best.get(col, 0.0), 6) if pd.notna(best.get(col)) else 0.0
        row['Distance to positive'] = round(best.get('TSum', 0.0), 6)
        row['Distance to negative'] = round(best.get('FSum', 0.0), 6)

        rows.append(row)

    if rows:
        summary_df = pd.DataFrame(rows)
        summary_df.to_csv(output_file, index=False)
        print(f"\nSummary saved to {output_file} ({len(rows)} ligands)")
    else:
        print("\n[WARNING] No matching results found for summary.")


# ============================================================
# PIPELINE LOGIC (COORDINATE MATH & CONFIG)
# ============================================================

def generate_random_complex(receptor_path, ligand_path, output_path):
    """
    Calculates centers, applies random translation, and writes a combined PDB.
    """
    # Read receptor coordinates
    receptor_coords = []
    with open(receptor_path, "r") as f:
        for line in f:
            if line.startswith(("ATOM","HETATM")):
                receptor_coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))

    cx = sum(c[0] for c in receptor_coords) / len(receptor_coords)
    cy = sum(c[1] for c in receptor_coords) / len(receptor_coords)
    cz = sum(c[2] for c in receptor_coords) / len(receptor_coords)
    
    print(f"Receptor center: ({cx:.2f}, {cy:.2f}, {cz:.2f})")

    # Read ligand coordinates
    ligand_coords = []
    with open(ligand_path, "r") as f:
        for line in f:
            if line.startswith(("ATOM","HETATM")):
                ligand_coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))

    lx = sum(c[0] for c in ligand_coords) / len(ligand_coords)
    ly = sum(c[1] for c in ligand_coords) / len(ligand_coords)
    lz = sum(c[2] for c in ligand_coords) / len(ligand_coords)
    
    print(f"Ligand center: ({lx:.2f}, {ly:.2f}, {lz:.2f})")

    # Compute translation
    distance = random.uniform(10, 15)
    dx = cx - lx + distance
    dy = cy - ly + random.uniform(-5, 5)
    dz = cz - lz + random.uniform(-5, 5)

    print("Translation vector:")
    print(f"dx = {dx:.3f}")
    print(f"dy = {dy:.3f}")
    print(f"dz = {dz:.3f}")

    # Translate ligand lines
    with open(ligand_path, "r") as f:
        ligand_lines = f.readlines()
    ligand_translated = []
    for line in ligand_lines:
        if line.startswith(("ATOM","HETATM")):
            x = float(line[30:38]) + dx
            y = float(line[38:46]) + dy
            z = float(line[46:54]) + dz
            new_line = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
            ligand_translated.append(new_line)
        else:
            ligand_translated.append(line)

    # Write complex PDB
    with open(output_path, "w") as out:
        with open(receptor_path, "r") as f:
            out.write(f.read())
        out.write("\n")
        out.writelines(ligand_translated)


def load_config():
    config = configparser.ConfigParser()
    config.read(BASE_DATA / "config")
    try:
        receptor_name = config["Input_Files"]["pdb_receptor"]

        # Discover region IDs and compute totals dynamically from [Regions] section
        region_ids = set()
        region_totals = {}
        if "Regions" in config:
            for _, val in config["Regions"].items():
                parts = [p.strip() for p in val.split(",")]
                if len(parts) == 3:
                    try:
                        rid = int(parts[0])
                        weight = float(parts[2])
                        region_ids.add(rid)
                        region_totals[rid] = region_totals.get(rid, 0.0) + weight
                    except ValueError:
                        pass

        # --- Hardcoded Rosetta option lists ---
        cst_weight = config.getfloat("CST_Weight", "cst_weight", fallback=0.1)
        global_opts = [
            "-docking", "-partners", "A_B",
            "-docking:randomize1", "-docking:randomize2", "-docking:spin",
            "-docking:use_ellipsoidal_randomization", "true",
            "-constraints:cst_file", "constraints.cst",
            "-constraints:cst_weight", str(cst_weight),
            "-ex1", "-ex2aro", "-out:pdb", "-overwrite",
        ]
        local_opts = [
            "-docking", "-partners", "A_B",
            "-docking:docking_local_refine", "true",
            "-ex1", "-ex2aro", "-out:file:fullatom", "-overwrite",
        ]

        return {
            "i_global_max": int(config["General_Constants"]["i_global_max"]),
            "i_local_max": int(config["General_Constants"]["i_local_max"]),
            "ratio": config.getfloat("Reward", "ratio", fallback=1.0),
            "ratio_increment": config.getfloat("Reward", "ratio_increment", fallback=0.0),
            "success_cycles": config.getint("General_Constants", "success_cycles", fallback=1),
            "receptor": BASE_DATA / receptor_name,
            "global_opts": global_opts,
            "local_opts": local_opts,
            "region_ids": sorted(region_ids),
            "region_totals": {rid: region_totals[rid] for rid in sorted(region_ids)},
        }
    except KeyError as e:
        print(f"[FATAL] Missing config key {e}")
        sys.exit(1)


def write_runtime_ratio_config(ratio_value, cycle_tag):
    """Create a temporary config file for this cycle with an updated ratio."""
    source_config = BASE_DATA / "config"
    runtime_config = BASE_DATA / f"runtime_config_{cycle_tag}.ini"

    parser = configparser.ConfigParser()
    parser.read(source_config)
    if "General_Constants" not in parser:
        parser["General_Constants"] = {}

    parser["General_Constants"]["ratio"] = f"{ratio_value:.10g}"

    with open(runtime_config, "w", encoding="utf-8") as f:
        parser.write(f)

    return runtime_config


# ============================================================
# STAGE WRAPPERS
# ============================================================

def run_preprocessing():
    print("\n===== PREPROCESSING =====")
    try:
        subprocess.run([sys.executable, SCRIPT_PRE], check=True)
    except subprocess.CalledProcessError:
        print("[WARNING] Preprocessing failed. Attempting to continue...")


def run_rosetta(stage, ligand_name, iteration, input_file, options, pdbs_dir, scores_dir, nstruct=1):
    prefix = f"{stage}_{ligand_name}_t{iteration}_"
    score = f"score_{stage}_{ligand_name}_t{iteration}.sc"
    pdbs_dir.mkdir(parents=True, exist_ok=True)
    scores_dir.mkdir(parents=True, exist_ok=True)
    prefix_path = pdbs_dir / prefix
    score_path = scores_dir / score

    cmd = [
        "docker","run","--rm","--platform","linux/amd64",
        "-v",f"{HOST_DATA}:/data",
        "-w","/data",
        ROSETTA_IMAGE,
        "bash","-lc",
        f"rosettadock {' '.join(options)} "
        f"-in:file:s {input_file} "
        f"-out:prefix {prefix_path} "
        f"-out:file:scorefile {score_path} "
        f"-nstruct {nstruct}"
    ]

    run_command(cmd, f"Rosetta {stage}")
    pdb_path = rename_rosetta_output(
        prefix,
        f"{stage}_{ligand_name}_t{iteration}.pdb",
        output_dir=pdbs_dir,
        search_dir=pdbs_dir,
    )

    return pdb_path, score_path


def run_pisa(ligand_name, iteration, pdb_file, ccp4_dir):
    print(f"Iteration {iteration}: PISAePDB")
    ccp4_dir.mkdir(parents=True, exist_ok=True)
    folder_in = ccp4_dir / f"ccp4_in_{ligand_name}_t{iteration}"
    folder_out = ccp4_dir / f"ccp4_out_{ligand_name}_t{iteration}"

    folder_in.mkdir(exist_ok=True)
    folder_out.mkdir(exist_ok=True)

    try:
        with open(pdb_file, "r") as f:
            filtered = [line for line in f if line.startswith("ATOM") or line.startswith("TER")]

        safe_pdb_name = build_pisa_input_name(pdb_file.name)
        #if safe_pdb_name != pdb_file.name:
            #print(f"[INFO] PISA input renamed to safe filename: {pdb_file.name} -> {safe_pdb_name}")

        pdb_input = folder_in / safe_pdb_name
        with open(pdb_input,"w") as f:
            f.writelines(filtered)

        cmd = [
            "docker","run","--rm",
            "-v",f"{HOST_DATA}:/data",
            "-w","/data",
            CCP4_IMAGE,
            "bash","-c",
            f"/ccp4/bin/run {folder_in.relative_to(BASE_DATA)} {folder_out.relative_to(BASE_DATA)}"
        ]
        run_command(cmd, "CCP4/PISA")
    except Exception as e:
        print(f"[ERROR] PISA helper failed: {e}")

    return folder_in, folder_out


# ============================================================
# MAIN PIPELINE
# ============================================================

def main():
    run_preprocessing() # RESTORED
    cfg = load_config()
    ensure_log(cfg["region_ids"])

    ligand_dir = BASE_DATA / "ligands"
    if not ligand_dir.exists():
        print("[FATAL] No ligand folder found.")
        sys.exit(1)

    ligands = list(ligand_dir.glob("*.pdb"))
    if not ligands:
        print("[FATAL] No ligand files found.")
        sys.exit(1)

    for ligand in ligands:
        name = ligand.stem
        print(f"\n=== PROCESSING {ligand.name} ===")

        # Create output structure for this ligand
        ligand_dir = BASE_DATA / "final_results" / "PDBs" / name
        pdbs_dir = ligand_dir / "useful_files"
        scores_dir = ligand_dir / "scores"
        ccp4_dir = ligand_dir / "ccp4"
        approved_dir = ligand_dir / "approved_complexes"

        for d in [pdbs_dir, scores_dir, ccp4_dir, approved_dir]:
            d.mkdir(parents=True, exist_ok=True)

        complex_file = pdbs_dir / f"c_{name}.pdb"
        try:
            generate_random_complex(cfg["receptor"], ligand, complex_file)
        except Exception as e:
            print(f"[ERROR] Could not generate complex for {name}: {e}")
            continue

        start_complex = complex_file
        global_iteration = 0
        max_cycles = max(1, cfg["success_cycles"])

        print("Region totals (sum of weights from config):")
        for rid, total in cfg["region_totals"].items():
            print(f"Region {rid}: {total:.6f}")

        for cycle in range(1, max_cycles + 1):
            ratio_value = cfg["ratio"] + (cycle - 1) * cfg["ratio_increment"]
            run_name = f"{name}_c{cycle}"
            
            # --- INTEGRATED SELECTION & CALCULATION LOGIC ---

            current_local_opts = cfg["local_opts"]

            print(f"\n--- Cycle {cycle}/{max_cycles} ---")
            print(f"Ratio: {ratio_value:.4g}")
            if cycle == 1:
                print(f"[CMD] Rosetta global opts: {' '.join(cfg['global_opts'])}")
            print(f"[CMD] Rosetta local opts:  {' '.join(current_local_opts)}")

            runtime_config = write_runtime_ratio_config(ratio_value, f"{run_name}_{int(time.time())}")
            cycle_success = False
            approved_complex = None

            try:
                # Iteration Loop for this cycle
                cycle_n_max = cfg["i_global_max"] if cycle == 1 else cfg["i_local_max"]
                for iteration in range(1, cycle_n_max + 1):
                    global_iteration += 1
                    start = time.time()
                    current = start_complex
                    docker_error = False

                    
                    if cycle == 1:
                        print(f"Iteration {iteration}: Rosetta global")
                        current, score = run_rosetta("global", run_name, global_iteration, start_complex, cfg["global_opts"], pdbs_dir, scores_dir)
                        if current is None:
                            docker_error = True
                    else:
                        current = start_complex

                    if not docker_error:
                        print(f"Iteration {iteration}: Rosetta local")
                        current, score = run_rosetta("local", run_name, global_iteration, current, current_local_opts, pdbs_dir, scores_dir)
                        if current is None:
                            docker_error = True

                    if docker_error:
                        log_result(start_complex.name, global_iteration, "Failed", round(time.time() - start, 2), region_ids=cfg["region_ids"])
                        continue

                    folder_in, folder_out = run_pisa(run_name, iteration, current, ccp4_dir)
                    elapsed = round(time.time() - start, 2)
                    success, data, tsv_file, tsv_files = evaluate_pisa_outputs(folder_out, runtime_config)

                    if not tsv_files:
                        log_result(current.name, global_iteration, "No Interaction", elapsed, region_ids=cfg["region_ids"])
                        shutil.rmtree(folder_in)
                        shutil.rmtree(folder_out)
                        continue

                    if data is None:
                        log_result(current.name, global_iteration, "No Valid Interaction Data", elapsed, region_ids=cfg["region_ids"])
                        shutil.rmtree(folder_in)
                        shutil.rmtree(folder_out)
                        continue

                    if success:
                        print("Structure approved! Saving solution and continuing cascade...")
                        approved_elapsed = round(time.time() - start, 2)
                        log_result(current.name, global_iteration, "Approved", approved_elapsed, data, cfg["region_ids"])

                        # Save approved complex
                        approved_copy = approved_dir / f"cycle{cycle}_iter{global_iteration}_{current.name}"
                        shutil.copy2(current, approved_copy)

                        # Extract structure summaries from interface TSV
                        if tsv_file is not None:
                            extract_structure_summaries(tsv_file, ligand_dir)

                        if folder_in.exists():
                            shutil.rmtree(folder_in)

                        cycle_success = True
                        approved_complex = current
                        break

                    print("Structure rejected")
                    log_result(current.name, global_iteration, "Rejected", elapsed, data, cfg["region_ids"])
                    shutil.rmtree(folder_in)
                    shutil.rmtree(folder_out)

                    # Iteration cleanup
                    pattern = f"*_t{global_iteration}.*"
                    for search_dir in [pdbs_dir, scores_dir]:
                        for f in search_dir.glob(pattern):
                            try:
                                f.unlink()
                            except Exception:
                                pass
            finally:
                if runtime_config.exists():
                    runtime_config.unlink()

            if not cycle_success:
                print(f"No approved structure found in cycle {cycle}. Stopping cascade for {ligand.name}.")
                break

            start_complex = approved_complex

        print(f"Finished {ligand.name}")

    # Generate summary from control results
    control_csv = ensure_control_csv(BASE_DATA)
    if control_csv:
        summary_file = BASE_DATA / "final_results" / "summary.csv"
        generate_summary(control_csv, LOG_FILE, summary_file, cfg["region_ids"])
    else:
        print("\n[INFO] No control_results file found. Skipping summary.")

        



if __name__ == "__main__":
    main()