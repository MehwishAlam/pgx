import pandas as pd
from pathlib import Path
import json

def filter_pgx_by_sample(sample_id: str) -> pd.DataFrame:
    """
    Load 'PGX OA Genotyping.xls' (QuantStudio export) and filter rows by Sample ID.

    Returns a DataFrame with columns:
        ['Gene Symbol', 'NCBI SNP Reference', 'Sample ID', 'Call']
    """
    pgx_xls_path = "data_input/08212025 PGX OA Genotyping Data (1) (1) (1).xlsx"
    pgx_path = Path(pgx_xls_path)

    # 1) Read entire sheet with no header (because there are metadata rows above the table)
    # For .xls, pandas normally uses the xlrd engine automatically.
    raw = pd.read_excel(pgx_path, sheet_name=0, header=None)

    # 2) Find the header row: first row whose first cell is "Assay Name"
    header_row_idx = None
    for i in range(len(raw)):
        if str(raw.iloc[i, 0]).strip() == "Assay Name":
            header_row_idx = i
            break

    if header_row_idx is None:
        raise ValueError("Could not find header row starting with 'Assay Name' in the XLS file.")

    # 3) Build the real dataframe: use that row as header
    header = raw.iloc[header_row_idx].tolist()
    df = raw.iloc[header_row_idx + 1 :].copy()  # data below header
    df.columns = header
    df.reset_index(drop=True, inplace=True)

    # 4) Keep only required columns
    cols = ["Gene Symbol", "NCBI SNP Reference", "Sample ID", "Call"]
    for c in cols:
        if c not in df.columns:
            raise KeyError(f"Required column '{c}' not found in XLS file.")

    df = df[cols]

    # 5) Filter by Sample ID (string-compare, stripped)
    sample_id_str = str(sample_id).strip()
    df = df[df["Sample ID"].astype(str).str.strip() == sample_id_str]

    # 6) Remove duplicate wells for the same SNP+Sample (if any)
    df = df.drop_duplicates(subset=["NCBI SNP Reference", "Sample ID"])

    return df

def get_unique_genes(filtered_df: pd.DataFrame) -> list:
    """
    Return a list of unique gene symbols from the filtered PGX dataframe.
    """
    return sorted(filtered_df["Gene Symbol"].drop_duplicates().tolist())

def summarize_genes_with_json(
    unique_genes: list,
    filtered_df: pd.DataFrame,
    json_folder: str = "KG/allele_definition/json_file",
):
    """
    For each gene in unique_genes:
      - Extract all rows (SNP + Call) from filtered_df
      - Print a small summary
      - Load the corresponding JSON file from json_folder/<GENE>.json

    Returns:
        dict[gene] = {
            "calls": [
                {"snp": ..., "call": ...},
                ...
            ],
            "allele_def": <loaded JSON dict or None>
        }
    """
    json_base = Path(json_folder)
    summary = {}

    for gene in unique_genes:
        gene_rows = filtered_df[filtered_df["Gene Symbol"] == gene]

        print(f"\n===== Gene: {gene} =====")
        gene_calls = []
        for _, row in gene_rows.iterrows():
            snp = row["NCBI SNP Reference"]
            call = row["Call"]
            gene_calls.append({"snp": snp, "call": call})
            print(f"SNP: {snp}  |  Call: {call}")

        summary[gene] = {
            "calls": gene_calls
        }

    return summary

def map_genotypes_to_star_alleles(
    gene_summary: dict,
    json_folder: str = "KG/allele_definition/json_file",
):
    """
    Custom mapping rule:
    - Find allele for second base first (e.g. 'G' in A/G)
    - Then take the closest matching allele *above it* for the first base
    """
    base_dir = Path(json_folder)
    all_results = {}

    for gene, info in gene_summary.items():
        json_path = base_dir / f"{gene}_allele_definition_table.json"
        if not json_path.exists():
            print(f"⚠ JSON missing for {gene}")
            continue

        with json_path.open("r", encoding="utf-8") as f:
            allele_def = json.load(f)

        rsids = allele_def.get("rsID")
        allele_keys = [k for k in allele_def.keys() if k.startswith("*")]

        gene_results = []

        for call_entry in info.get("calls", []):
            snp = call_entry["snp"]
            call = call_entry["call"]

            # skip failed calls
            if call.upper() in {"NOAMP", "NOCALL", "UND", "INVALID"}:
                continue

            if snp not in rsids:
                continue

            pos = rsids.index(snp)
            a1, a2 = call.replace(" ", "").split("/")

            # Step 1 → find second allele match first (a2)
            second_star = None
            second_idx = None
            for idx, star in enumerate(allele_keys):
                values = allele_def[star]
                base = str(values[pos])
                if base == a2:
                    second_star = star
                    second_idx = idx
                    break

            if second_star is None:
                continue

            # Step 2 → scan backwards to find nearest match for a1
            first_star = None
            for idx in range(second_idx - 1, -1, -1):
                star = allele_keys[idx]
                base = str(allele_def[star][pos])
                if base == a1:
                    first_star = star
                    break

            # If no match above, fallback to any match in list
            if first_star is None:
                for idx, star in enumerate(allele_keys):
                    base = str(allele_def[star][pos])
                    if base == a1:
                        first_star = star
                        break

            diplotype = f"{first_star}/{second_star}"

            result = {
                "snp": snp,
                "call": call,
                "position_index": pos,
                "first_star": first_star,
                "second_star": second_star,
                "diplotype": diplotype
            }
            gene_results.append(result)

            print(f"[{gene}] {snp} {call} → {diplotype}")

        all_results[gene] = gene_results

    return all_results

def map_star_functions_from_json(
    star_alleles: dict,
    json_folder: str = "KG/allele_functionality/json_file"
):
    """
    Load star allele functionality from per-gene JSON files:
      KG/allele_functionality/json_file/<GENE>.json

    Example JSON format:
      {
        "CYP3A5": {
          "*1": {
            "Allele Clinical Functional Status (Required)": "Normal function",
            ...
          },
          "*3": {
            "Allele Clinical Functional Status (Required)": "No function",
            ...
          }
        }
      }

    Returns:
      {
        "CYP3A5": {
          "*1": "Normal function",
          "*3": "No function"
        }
      }
    """
    base_dir = Path(json_folder)
    final_result = {}

    for gene, entries in star_alleles.items():
        json_path = base_dir / f"{gene}.json"

        if not json_path.exists():
            print(f"⚠ Functionality JSON not found for {gene}: {json_path}")
            final_result[gene] = {}
            continue

        # Load gene functionality JSON
        with open(json_path, "r", encoding="utf-8") as f:
            gene_json = json.load(f)

        gene_map = gene_json.get(gene, {})
        gene_result = {}

        # Extract star alleles found in patient's result
        stars_present = set()
        for entry in entries:
            if entry.get("first_star"):
                stars_present.add(entry["first_star"])
            if entry.get("second_star"):
                stars_present.add(entry["second_star"])

        for star in sorted(stars_present):
            func_info = gene_map.get(star)

            if isinstance(func_info, dict):
                func = func_info.get(
                    "Allele Clinical Functional Status (Required)",
                    "Unknown function"
                )
                gene_result[star] = func
            else:
                gene_result[star] = "Unknown function"

        final_result[gene] = gene_result

    return final_result

def map_diplotypes_to_phenotypes(
    star_alleles: dict,
    star_functions: dict,
    phenotype_json_folder: str = "KG/diplotype-phenotype/json_file",
):
    """
    For each gene in star_alleles:
      - collect unique diplotypes from entries (entry["diplotype"])
      - load <GENE>_diplotype_phenotype.json
      - map each diplotype -> Phenotype, Activity Score, EHR Priority
      - attach star-level functions for the two stars

    Parameters
    ----------
    star_alleles : dict
        e.g. {
          "CYP3A5": [
            {"snp": "...", "call": "...", "first_star": "*1", "second_star": "*3", "diplotype": "*1/*3"},
            {"snp": "...", "call": "...", "first_star": "*1", "second_star": "*1", "diplotype": "*1/*1"},
            ...
          ],
          "CYP2C19": [...],
          ...
        }

    star_functions : dict
        e.g. {
          "CYP3A5": {
            "*1": "Normal function",
            "*3": "No function",
            "*6": "No function",
            ...
          },
          "CYP2C19": {...}
        }

    phenotype_json_folder : str
        Folder where gene-level phenotype JSON files live:
        KG/diplotype_phenotype/json_file/<GENE>_diplotype_phenotype.json

    Returns
    -------
    dict
        {
          "CYP3A5": {
            "*1/*1": {
              "alleles": ["*1", "*1"],
              "allele_functions": {"*1": "Normal function"},
              "phenotype": "CYP3A5 Normal Metabolizer",
              "activity_score": "",
              "ehr_priority": "Abnormal/Priority/High Risk",
              "source_file": "CYP3A5_diplotype_phenotype.json"
            },
            "*1/*3": {
              "alleles": ["*1", "*3"],
              "allele_functions": {
                "*1": "Normal function",
                "*3": "No function"
              },
              "phenotype": "CYP3A5 Intermediate Metabolizer",
              "activity_score": "",
              "ehr_priority": "Abnormal/Priority/High Risk",
              "source_file": "CYP3A5_diplotype_phenotype.json"
            },
            ...
          },
          "CYP2C19": {
            ...
          }
        }
    """
    base_dir = Path(phenotype_json_folder)
    result = {}

    for gene, entries in star_alleles.items():
        # --- 1) collect unique diplotypes for this gene ---
        diplotypes = set()
        for e in entries:
            d = e.get("diplotype")
            if not d:
                # fall back to first/second_star if needed
                s1, s2 = e.get("first_star"), e.get("second_star")
                if s1 and s2:
                    d = f"{s1}/{s2}"
                else:
                    continue
            diplotypes.add(str(d).strip())

        if not diplotypes:
            result[gene] = {}
            continue

        # --- 2) load phenotype JSON for this gene ---
        phen_file = base_dir / f"{gene}_diplotype_phenotype.json"
        gene_phens = {}
        phen_json = {}
        if phen_file.exists():
            with phen_file.open("r", encoding="utf-8") as f:
                phen_json = json.load(f)
            gene_block = phen_json.get(gene, {})
        else:
            print(f"⚠ Phenotype JSON not found for {gene}: {phen_file}")
            gene_block = {}

        gene_func_map = star_functions.get(gene, {})

        # --- 3) map each diplotype to phenotype info ---
        for d in sorted(diplotypes):
            # parse stars from diplotype string
            parts = d.replace(" ", "").split("/")
            if len(parts) != 2:
                continue
            s1, s2 = parts[0], parts[1]

            # lookup phenotype in JSON
            phen_info = gene_block.get(d)
            if phen_info is None and gene_block:
                # try reversed
                rev = f"{s2}/{s1}"
                phen_info = gene_block.get(rev)

            if phen_info is None:
                phenotype = "Unknown"
                activity_score = ""
                ehr_priority = ""
            else:
                phenotype = phen_info.get("Phenotype", "Unknown")
                activity_score = phen_info.get("Activity Score", "")
                ehr_priority = phen_info.get("EHR Priority Notation", "")

            # allele functions
            allele_functions = {}
            for star in {s1, s2}:
                allele_functions[star] = gene_func_map.get(star, "Unknown function")

            gene_phens[d] = {
                "alleles": [s1, s2],
                "allele_functions": allele_functions,
                "phenotype": phenotype,
                "activity_score": activity_score,
                "ehr_priority": ehr_priority,
                "source_file": phen_file.name if phen_file.exists() else None,
            }

        result[gene] = gene_phens

    return result

sample = "EDX2508083837"

def run_pgx_technical_report(sample: str, return_all_steps: bool = False):
    result = filter_pgx_by_sample(sample)

    unique_gene_list = get_unique_genes(result)
    print(" Unique Gene List : ",unique_gene_list)

    gene_summary = summarize_genes_with_json(unique_gene_list, result)
    print(" Gene Summary : ",gene_summary)

    star_alleles = map_genotypes_to_star_alleles(gene_summary)
    print(" Star Alleles : ",star_alleles)

    star_functions = map_star_functions_from_json(star_alleles)
    print(" Star Functions : ",star_functions)

    gene_phenotypes = map_diplotypes_to_phenotypes(star_alleles, star_functions)
    print(" Gene Phenotypes : ",gene_phenotypes)
    
    if return_all_steps:
        return {
            "filtered_data": result,
            "unique_gene_list": unique_gene_list,
            "gene_summary": gene_summary,
            "star_alleles": star_alleles,
            "star_functions": star_functions,
            "gene_phenotypes": gene_phenotypes
        }
    
    return gene_phenotypes