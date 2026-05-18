"""
scRNA-seq analysis for TRPM5-KLF4 hypothesis
Uses CELLxGENE API to query human lung scRNA-seq data.

Tests if tuft cell markers and basal cell KLF4 show cell-type-specific expression.
"""

import requests
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import time
import json

API = "https://api.cellxgene.cziscience.com"

# Genes of interest
TUFT_GENES = ["POU2F3", "TRPM5", "GNAT3", "CHAT", "IL25", "AVIL"]
BASAL_GENES = ["KRT5", "TP63", "KRT14", "NGFR"]
CHANNEL_RECEPTORS = ["CHRM3", "CHRNA3", "CHRNA5", "CHRNA7", "IL17RB", "CYSLTR1"]
DOWNSTREAM = ["KLF4", "KRT13", "HIF1A"]

ALL_GENES = list(dict.fromkeys(TUFT_GENES + BASAL_GENES + CHANNEL_RECEPTORS + DOWNSTREAM))


def find_best_lung_dataset():
    """Find the best human lung dataset from CELLxGENE."""
    print("[1] Fetching all CELLxGENE collections...")
    r = requests.get(f"{API}/curation/v1/collections", timeout=60)
    collections = r.json()
    print(f"    Total collections: {len(collections)}")

    lung_ds = []
    for col in collections:
        for ds in col.get("datasets", []):
            # Check human
            organisms = [o.get("label", "") for o in ds.get("organism", [])]
            if "Homo sapiens" not in organisms:
                continue

            # Check lung tissue
            tissues = [t.get("label", "").lower() for t in ds.get("tissue", [])]
            tissue_str = " ".join(tissues)
            lung_keywords = ["lung", "respiratory", "airway", "bronch", "trachea",
                           "nasal epithelium", "nasal mucosa"]
            is_lung = any(kw in tissue_str for kw in lung_keywords)

            if not is_lung:
                continue

            cell_count = ds.get("cell_count", 0)
            if cell_count < 1000:
                continue

            disease = [d.get("label", "").lower() for d in ds.get("disease", [])]
            is_normal = "normal" in " ".join(disease)

            assay = [a.get("label", "") for a in ds.get("assay", [])]
            assay_str = " ".join(assay)

            lung_ds.append({
                "collection_name": col.get("name", "")[:120],
                "dataset_id": ds.get("dataset_id", ""),
                "label": ds.get("label", "")[:120],
                "cell_count": cell_count,
                "tissue": tissue_str[:80],
                "disease": " ".join(disease)[:60],
                "assay": assay_str[:60],
                "is_normal": is_normal,
                "collection_id": col.get("collection_id", ""),
                "doi": col.get("doi", ""),
            })

    # Prefer large, normal lung datasets
    lung_ds.sort(key=lambda x: (x["is_normal"], x["cell_count"]), reverse=True)
    return lung_ds


def get_dataset_wms(dataset_id):
    """Fetch WMS dataset metadata including cell types."""
    try:
        r = requests.get(f"{API}/wms/v2/dataset/{dataset_id}", timeout=30)
        if r.status_code == 200:
            return r.json()
    except Exception:
        pass
    return None


def query_expression_wms(dataset_id, gene_symbol):
    """Query gene expression summary from WMS API.

    Returns list of {cell_type: str, expression: float} or None.
    """
    url = f"{API}/wms/v2/expression"
    payload = {
        "dataset_id": dataset_id,
        "gene_symbols": [gene_symbol],
    }
    try:
        r = requests.post(url, json=payload, timeout=60)
        if r.status_code == 200:
            result = r.json()
            return parse_expression_response(result, gene_symbol)
        elif r.status_code == 404:
            return None
        else:
            # Try GET-based endpoint
            alt_url = f"{API}/wms/v2/expression/summary/{gene_symbol}"
            alt_params = {"dataset_id": dataset_id}
            r2 = requests.get(alt_url, params=alt_params, timeout=60)
            if r2.status_code == 200:
                return parse_expression_response(r2.json(), gene_symbol)
            return None
    except Exception as e:
        return None


def parse_expression_response(data, gene_symbol):
    """Parse WMS expression response into a simple format."""
    # The response format varies; try common patterns
    if isinstance(data, dict):
        if "expression" in data:
            return {"gene": gene_symbol, "data": data["expression"]}
        if "results" in data:
            return {"gene": gene_symbol, "data": data["results"]}
        return {"gene": gene_symbol, "data": data}
    return {"gene": gene_symbol, "data": data}


# ============================================================
# Fallback: Use CELLxGENE Explorer's CSV export
# ============================================================
EXPLORER_URL = "https://cellxgene.cziscience.com"


def get_explorer_gene_data(dataset_id, gene):
    """Try the Explorer's gene expression endpoint."""
    url = f"{EXPLORER_URL}/api/v0.2/dataset/{dataset_id}/expression/gene/{gene}"
    try:
        r = requests.get(url, timeout=30)
        if r.status_code == 200:
            return r.json()
    except Exception:
        pass
    return None


# ============================================================
# Alternative: Direct download from known lung atlas
# ============================================================

# The HLCA core dataset is available at this known URL.
# We can use the processed data from the CELLxGENE data portal.
# For now, let's focus on API-based queries.

HLCA_DATASET_ID = "8aa49acf-4e7d-4747-a3f2-a15a421f020d"  # HLCA main dataset
HLCA_COLLECTION_ID = "f5e2c0e0-1b6c-4b2a-bf6b-b0208b1e4c96"


def try_get_expression_via_portal(dataset_id, gene_list):
    """Try the CELLxGENE data portal expression per cell type.

    Uses the known pattern from CELLxGENE Explorer.
    """
    results = {}
    for gene in gene_list:
        # Explorer API endpoint for gene expression
        url = f"{EXPLORER_URL}/d/{dataset_id}/api/v0.3/expression/gene/{gene}"
        try:
            r = requests.get(url, timeout=30)
            if r.status_code == 200:
                results[gene] = r.json()
            else:
                # Alternative endpoint
                url2 = f"https://api.cellxgene.cziscience.com/curation/v1/datasets/{dataset_id}/expression/{gene}"
                r2 = requests.get(url2, timeout=30)
                if r2.status_code == 200:
                    results[gene] = r2.json()
        except Exception as e:
            pass
        time.sleep(0.2)
    return results


# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    print("=" * 65)
    print("scRNA-seq: Tuft Cell -> Basal Cell Communication")
    print("=" * 65)

    # Step 1: Find lung datasets
    lung_datasets = find_best_lung_dataset()
    print(f"\n  Found {len(lung_datasets)} human lung scRNA-seq datasets.\n")

    for ds in lung_datasets[:10]:
        status = "NORMAL" if ds["is_normal"] else ds["disease"]
        print(f"  [{ds['cell_count']:>8,} cells] {ds['collection_name'][:70]}")
        print(f"    ID: {ds['dataset_id']}")
        print(f"    Tissue: {ds['tissue']} | {status}")
        print(f"    Assay: {ds['assay'][:50]}")
        print()

    if not lung_datasets:
        print("  ERROR: No lung datasets found!")
        exit(1)

    # Step 2: Pick best dataset
    best = lung_datasets[0]
    dataset_id = best["dataset_id"]
    print(f"[2] Selected: {best['collection_name'][:80]}")
    print(f"    Cells: {best['cell_count']:,} | Tissue: {best['tissue']}")
    print(f"    Dataset ID: {dataset_id}")

    # Step 3: Get WMS metadata (cell types)
    print(f"\n[3] Getting dataset metadata...")
    wms = get_dataset_wms(dataset_id)
    if wms:
        print(f"    Title: {wms.get('title', '?')}")
        schema = wms.get("schema", {})
        annotations = schema.get("annotations", {})
        for ann in annotations.get("obs", []):
            name = ann.get("name", "")
            if name in ("cell_type", "cell_type_ontology_term_id",
                       "celltype", "CellType", "Cell Type"):
                categories = ann.get("categories", [])
                print(f"    Cell types ({len(categories)}):")
                for ct in sorted(categories)[:30]:
                    print(f"      - {ct}")
    else:
        print("    Could not get WMS metadata.")

    # Step 4: Try querying expression
    print(f"\n[4] Trying expression queries...")

    # Try all known endpoints for a few key genes
    for gene in ["TRPM5", "POU2F3", "KLF4"]:
        print(f"  {gene}:")
        # WMS POST
        result = query_expression_wms(dataset_id, gene)
        if result:
            print(f"    WMS POST: {str(result)[:200]}")
        else:
            print(f"    WMS POST: No result")

        # Explorer API
        exp_data = get_explorer_gene_data(dataset_id, gene)
        if exp_data:
            print(f"    Explorer: {str(exp_data)[:200]}")
        else:
            print(f"    Explorer: No result")

        time.sleep(0.3)

    # Step 5: Try with known HLCA dataset ID as fallback
    print(f"\n[5] Trying HLCA dataset directly ({HLCA_DATASET_ID})...")
    for gene in ["TRPM5", "POU2F3", "KLF4"]:
        result = query_expression_wms(HLCA_DATASET_ID, gene)
        if result:
            print(f"  {gene}: {str(result)[:300]}")

    print("\n" + "=" * 65)
    print("Query complete.")
    print("=" * 65)
