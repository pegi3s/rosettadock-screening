import csv
import configparser


def _load_config(config_file_path):
    config = configparser.ConfigParser()
    files_read = config.read(config_file_path)
    if not files_read:
        return None
    return config

# =========================================================
# 1. DICTIONARY EXTRACTION FUNCTION
# =========================================================
def extract_interface_dictionary(file_path):
    results_dict = {}
    try:
        row_8 = None
        row_11 = None
        with open(file_path, 'r', encoding='utf-8', newline='') as f:
            reader = csv.reader(f, delimiter='\t')
            for row_number, row in enumerate(reader, start=1):
                if row_number == 8:
                    row_8 = row
                elif row_number == 11:
                    row_11 = row
                    break

        if row_8 is None or row_11 is None:
            return {}

        col_limit = min(len(row_8), len(row_11))

        for i in range(col_limit):
                value_str = row_11[i].strip()
                key_str = row_8[i].strip()
                if not value_str or not key_str:
                    continue
                try:
                    value_num = float(value_str)
                    if value_num == 0.0 or value_num == 1.0:
                        try:
                            key = int(float(key_str))
                        except ValueError:
                            key = key_str
                        results_dict[key] = int(value_num)
                except ValueError:
                    continue

        return results_dict
    except FileNotFoundError:
        return {}
    except (OSError, csv.Error):
        return {}


def _get_region_ids(config):
    """Return sorted list of unique region IDs from [Regions] section."""
    if 'Regions' not in config:
        return []
    ids = set()
    for _, val in config['Regions'].items():
        parts = [p.strip() for p in val.split(',')]
        if len(parts) == 3:
            try:
                ids.add(int(parts[0]))
            except ValueError:
                pass
    return sorted(ids)


def _compute_region_totals(config):
    """Compute total weight sums per region from [Regions] section."""
    if 'Regions' not in config:
        return {}
    totals = {}
    for _, val in config['Regions'].items():
        parts = [p.strip() for p in val.split(',')]
        if len(parts) != 3:
            continue
        try:
            region = int(parts[0])
            weight = float(parts[2])
        except ValueError:
            continue
        totals[region] = totals.get(region, 0.0) + weight
    return totals


def _compute_region_sums(config, results_dict):
    """Returns dict {region_id: weighted_sum}."""
    if 'Regions' not in config:
        return {}
    sums = {}
    for _, value_string in config['Regions'].items():
        parts = [part.strip() for part in value_string.split(',')]
        if len(parts) != 3:
            continue
        try:
            region = int(parts[0])
            site = int(parts[1])
            weight = float(parts[2])
        except ValueError:
            continue
        site_found = site if site in results_dict else (str(site) if str(site) in results_dict else None)
        if site_found is not None:
            sums[region] = sums.get(region, 0.0) + weight * results_dict[site_found]
    return sums


def _compute_region_percentages(config, sums):
    """Returns dict {region_id: percentage}. Totals computed from [Regions] weights."""
    totals = _compute_region_totals(config)
    percs = {}
    for region, s in sums.items():
        total = totals.get(region, 0.0)
        percs[region] = (s / total * 100) if total != 0 else 0.0
    return percs


def _compute_tf_distances(config, percs):
    """Returns dict {region_id: (dist_t, dist_f)}. Reads T1{r} and F1{r} from [Variables_TF]."""
    if 'Variables_TF' not in config:
        return {region: (0.0, 0.0) for region in percs}
    dists = {}
    for region, p in percs.items():
        try:
            t = float(config['Variables_TF'][f'T1{region}'])
            f = float(config['Variables_TF'][f'F1{region}'])
            dists[region] = ((p - t) ** 2, (p - f) ** 2)
        except (KeyError, ValueError):
            dists[region] = (0.0, 0.0)
    return dists


def _compute_final_metric(config, tsum, fsum):
    try:
        ratio = float(config['General_Constants']['ratio'])
    except (KeyError, ValueError):
        ratio = 1.0

    if fsum == 0 or ratio == 0:
        return 0.0, False

    metric = tsum / (fsum / ratio)
    verdict = metric < 1
    return metric, verdict

# =========================================================
# 2. FUNCTION TO CALCULATE SUMS
# =========================================================
def calculate_region_sums(config_file_path, results_dict):
    """Returns dict {region_id: weighted_sum}."""
    config = _load_config(config_file_path)
    if config is None:
        return {}
    return _compute_region_sums(config, results_dict)

# =========================================================
# 3. FUNCTION TO CALCULATE PERCENTAGES
# =========================================================
def calculate_region_percentages(config_file_path, sums):
    """Returns dict {region_id: percentage}."""
    config = _load_config(config_file_path)
    if config is None:
        return {}
    return _compute_region_percentages(config, sums)

# =========================================================
# 4. FUNCTION TO CALCULATE DISTANCES (T and F)
# =========================================================
def calculate_tf_distances(config_file_path, percs):
    """Returns dict {region_id: (dist_t, dist_f)}."""
    config = _load_config(config_file_path)
    if config is None:
        return {}
    return _compute_tf_distances(config, percs)

# =========================================================
# 5. FUNCTION TO CALCULATE TOTAL SUMS (T and F)
# =========================================================
def calculate_total_sums(dists):
    """dists: dict {region_id: (dist_t, dist_f)}. Returns (tsum, fsum)."""
    tsum = sum(d[0] for d in dists.values())
    fsum = sum(d[1] for d in dists.values())
    return tsum, fsum

# =========================================================
# 6. FUNCTION TO CALCULATE THE METRIC AND FINAL RULE
# =========================================================
def calculate_final_metric(config_file_path, tsum, fsum):
    config = _load_config(config_file_path)
    if config is None:
        return 0.0, False
    return _compute_final_metric(config, tsum, fsum)

# =========================================================
# THE MAESTRO (FUNCTION CALLED BY THE MAIN PIPELINE)
# =========================================================
def evaluate(table_path, config_path="/data/config"):
    results_dict = extract_interface_dictionary(table_path)
    
    # IF TABLE READING FAILS, IT MUST RETURN 2 VALUES AT ONCE: False and None
    if not results_dict: 
        return False, None

    config = _load_config(config_path)
    if config is None:
        return False, None
        
    sums = _compute_region_sums(config, results_dict)
    percs = _compute_region_percentages(config, sums)
    dists = _compute_tf_distances(config, percs)
    tsum, fsum = calculate_total_sums(dists)

    metric, verdict = _compute_final_metric(config, tsum, fsum)

    # WE CREATE THE PACKAGE FOR EXCEL HERE (dynamic regions):
    math_data = {"Metrica": metric, "TSum": tsum, "FSum": fsum}
    for region_id in sorted(percs.keys()):
        math_data[f"P{region_id}"] = percs[region_id]
        math_data[f"S{region_id}"] = sums.get(region_id, 0.0)

    # WE ALWAYS RETURN 2 VALUES TO THE PIPELINE: The verdict and the data
    return verdict, math_data