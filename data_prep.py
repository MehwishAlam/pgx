import json
from pathlib import Path
import pandas as pd

def build_allele_jsons(
    excel_folder: str = "KG/allele_definition",
    json_folder: str = "KG/allele_definition/json_file",
) -> None:
    excel_dir = Path(excel_folder)
    out_dir = Path(json_folder)
    out_dir.mkdir(parents=True, exist_ok=True)

    for xls_path in excel_dir.iterdir():
        if xls_path.suffix.lower() not in {".xls", ".xlsx"}:
            continue

        df = pd.read_excel(xls_path, sheet_name=0, header=None)
        allele_dict = {}

        # ---------- find rsID row ----------
        rsid_row_idx = None
        for i in range(len(df)):
            first_col = str(df.iloc[i, 0]).strip()
            if first_col.lower() == "rsid":
                rsid_row_idx = i
                break

        if rsid_row_idx is None:
            print(f"⚠️ No rsID row in {xls_path.name} — skipping file")
            continue

        # ---------- find reference row (first usable row after rsID) ----------
        ref_row_idx = None
        for j in range(rsid_row_idx + 1, len(df)):
            key = str(df.iloc[j, 0]).strip()
            if not key:
                continue
            # skip the "<GENE> Allele" label row
            if key.lower().endswith("allele"):
                continue
            ref_row_idx = j
            break

        ref_values = df.iloc[ref_row_idx, 1:].tolist() if ref_row_idx is not None else None

        # ---------- build dictionary ----------
        for i in range(len(df)):
            key_raw = df.iloc[i, 0]
            if pd.isna(key_raw):
                continue

            key = str(key_raw).strip()
            if not key:
                continue

            # skip gene title row like "Gene:CYP2C9" or "Gene:VKORC1"
            if key.lower().startswith("gene:"):
                continue

            # skip "<GENE> Allele" label row
            if key.lower().endswith("allele"):
                continue

            # ----- rsID row -----
            if key.lower() == "rsid":
                vals = df.iloc[i, 1:].tolist()
                vals = [None if pd.isna(v) or str(v).strip() == "" else v for v in vals]
                allele_dict[key] = vals
                continue

            row_vals = df.iloc[i, 1:].tolist()

            # ----- rows AFTER rsID (genotype / allele rows) -----
            if ref_values is not None and i >= ref_row_idx:
                filled_vals = []
                for col_idx, v in enumerate(row_vals):
                    if pd.isna(v) or str(v).strip() == "":
                        v = ref_values[col_idx]  # copy from reference row
                    filled_vals.append(None if pd.isna(v) or str(v).strip() == "" else v)
                allele_dict[key] = filled_vals
                continue

            # ----- rows ABOVE rsID (properties like Common Name, Effect on protein, etc.) -----
            vals = [None if pd.isna(v) or str(v).strip() == "" else v for v in row_vals]
            allele_dict[key] = vals

        # ---------- save JSON ----------
        out_path = out_dir / f"{xls_path.stem}.json"
        with out_path.open("w", encoding="utf-8") as f:
            json.dump(allele_dict, f, ensure_ascii=False, indent=2)

        print(f"✔ JSON saved: {out_path}")

import pandas as pd
import json
from pathlib import Path


def build_functionality_json(
    ref_folder="KG/allele_functionality",
    json_output_folder="KG/allele_functionality/json_file"
):
    """
    Reads each *_allele_functionality_reference.xls file in ref_folder
    and saves as JSON with **all functionality columns preserved**
    """

    ref_dir = Path(ref_folder)
    out_dir = Path(json_output_folder)
    out_dir.mkdir(parents=True, exist_ok=True)

    for file in ref_dir.glob("*_allele_functionality_reference.xlsx"):

        print(f"📄 Processing: {file.name}")

        raw_df = pd.read_excel(file, sheet_name=0, header=None)

        # Detect GENE: row
        gene_name = None
        header_row = None

        for i in range(len(raw_df)):
            cell = str(raw_df.iloc[i, 0]).strip()
            if cell.startswith("GENE:"):
                gene_name = cell.split("GENE:")[-1].strip()
                header_row = i + 1
                break

        if gene_name is None:
            print(f"⚠ Skipped — GENE not found in {file.name}")
            continue

        # Load the actual table below header
        df = pd.read_excel(file, sheet_name=0, header=header_row)

        # Identify the allele column (first column)
        allele_col = df.columns[0]

        # Filter only star rows
        df = df[df[allele_col].astype(str).str.startswith("*")]

        # Build JSON structure:
        # GENE → { "*1": {col1:val, col2:val, ...}, "*3": {...}}
        gene_dict = {}

        for _, row in df.iterrows():
            star = str(row[allele_col]).strip()
            gene_dict[star] = {}

            for col in df.columns[1:]:
                gene_dict[star][col] = (
                    "" if pd.isna(row[col]) else str(row[col]).strip()
                )

        json_path = out_dir / f"{gene_name}.json"
        with open(json_path, "w", encoding="utf-8") as f:
            json.dump({gene_name: gene_dict}, f, indent=2)

        print(f"✔ JSON saved → {json_path}")

    print("\n🎯 All gene functionality files converted successfully!")

def convert_diplotype_phenotype_to_json(
    input_folder: str = "KG/diplotype-phenotype",
    output_folder: str = "KG/diplotype-phenotype/json_file"
) -> None:
    """
    Convert diplotype-phenotype Excel files in `input_folder` to JSON files in `output_folder`.

    Expected columns in each file:
      - "<GENE> Diplotype" (e.g., "CYP3A5 Diplotype")
      - "Activity Score"
      - "Coded Diplotype/Phenotype Summary"
      - "EHR Priority Notation"

    Output JSON format per gene file:
      {
        "<GENE>": {
          "<diplotype>": {
            "Activity Score": "...",
            "Phenotype": "...",
            "EHR Priority Notation": "..."
          },
          ...
        }
      }
    """
    in_dir = Path(input_folder)
    out_dir = Path(output_folder)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Process all .xls / .xlsx files in the folder
    for path in in_dir.glob("*.xls*"):
        print(f"📄 Processing diplotype file: {path.name}")

        df = pd.read_excel(path)

        if df.empty:
            print(f"⚠ File {path.name} is empty, skipping.")
            continue

        # First column header is like "CYP3A5 Diplotype"
        diplotype_col = df.columns[0]
        gene = str(diplotype_col).split()[0].strip()

        # Expected other columns
        activity_col = "Activity Score"
        summary_col = "Coded Diplotype/Phenotype Summary"
        ehr_col = "EHR Priority Notation"

        # Check required columns
        missing = [c for c in [activity_col, summary_col, ehr_col] if c not in df.columns]
        if missing:
            print(f"⚠ Missing columns {missing} in {path.name}, skipping.")
            continue

        gene_dict = {}

        for _, row in df.iterrows():
            diplotype = str(row[diplotype_col]).strip()
            if diplotype == "" or pd.isna(diplotype):
                continue

            activity_score = "" if pd.isna(row[activity_col]) else str(row[activity_col]).strip()
            phenotype = "" if pd.isna(row[summary_col]) else str(row[summary_col]).strip()
            ehr_priority = "" if pd.isna(row[ehr_col]) else str(row[ehr_col]).strip()

            gene_dict[diplotype] = {
                "Activity Score": activity_score,
                "Phenotype": phenotype,
                "EHR Priority Notation": ehr_priority
            }

        # Wrap in top-level gene key
        json_obj = {gene: gene_dict}

        out_path = out_dir / f"{gene}_diplotype_phenotype.json"
        with out_path.open("w", encoding="utf-8") as f:
            json.dump(json_obj, f, indent=2)

        print(f"✔ Saved JSON → {out_path}")

    print("\n🎯 Diplotype–phenotype conversion completed.")

# build_allele_jsons()
# build_functionality_json()
convert_diplotype_phenotype_to_json()