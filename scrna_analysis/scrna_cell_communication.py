"""
scRNA-seq Cell Communication Analysis for TRPM5-KLF4 Hypothesis
Uses GSE193816 (airway epithelial cells from endobronchial brushings).
Dataset: Allergic asthma vs control, in vivo human airway epithelium.

Key analysis:
  1. Identify brush/tuft cells and basal cells
  2. Validate marker gene expression patterns
  3. Run LIANA cell-cell communication inference
  4. Test: do tuft cells signal to basal cells via ACh-nAChR?
"""

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')

print("=" * 70)
print("scRNA-seq Cell Communication: Tuft Cell -> Basal Cell")
print("=" * 70)

# ============================================================
# 1. Load Data
# ============================================================
adata = sc.read_h5ad("D:/ai agent/first-cc/scrna_analysis/data/GSE193816_AEC_data.h5ad.gz")
print(f"\n[1] Data loaded: {adata.n_obs} cells x {adata.n_vars} genes")

# ============================================================
# 2. Explore metadata
# ============================================================
print(f"\n[2] Dataset metadata:")
print(f"    obs columns: {list(adata.obs.columns)}")
print(f"    var columns: {list(adata.var.columns)[:10]}")

# Check cell type annotations
for col in adata.obs.columns:
    if any(kw in col.lower() for kw in ['cell', 'type', 'cluster', 'annot', 'lineage']):
        print(f"    {col}: {sorted(adata.obs[col].unique().astype(str))[:30]}")
        print(f"    n_categories: {adata.obs[col].nunique()}")

# Check condition/group annotations
for col in adata.obs.columns:
    if any(kw in col.lower() for kw in ['condition', 'group', 'disease', 'treat', 'diagnosis', 'status', 'sample']):
        print(f"    {col}: {sorted(adata.obs[col].unique().astype(str))[:20]}")

# ============================================================
# 3. Genes of interest
# ============================================================
TUFT_MARKERS = ["TRPM5", "POU2F3", "GNAT3", "CHAT", "SLC18A3", "AVIL", "IL25", "LTC4S"]
BASAL_MARKERS = ["KRT5", "TP63", "KRT14", "NGFR", "EGFR"]
CHOLINERGIC_RECEPTORS = ["CHRM3", "CHRNA3", "CHRNA5", "CHRNA7", "CHRNB2", "CHRNB4"]
SQUAMOUS_MARKERS = ["KLF4", "KRT13", "HIF1A", "SOX2", "SPRR1A", "SPRR1B", "IVL"]
DOWNSTREAM = ["KRT13", "HIF1A", "TP63", "SOX2"]

ALL_GENES = list(dict.fromkeys(
    TUFT_MARKERS + BASAL_MARKERS + CHOLINERGIC_RECEPTORS + SQUAMOUS_MARKERS + DOWNSTREAM
))

# Check which genes are in the dataset
available_genes = [g for g in ALL_GENES if g in adata.var_names]
missing_genes = [g for g in ALL_GENES if g not in adata.var_names]
print(f"\n[3] Gene availability: {len(available_genes)}/{len(ALL_GENES)}")
if missing_genes:
    print(f"    Missing: {missing_genes}")
print(f"    Available: {available_genes}")

# ============================================================
# 4. Identify cell types via marker genes
# ============================================================
print(f"\n[4] Cell type identification via marker expression...")

# Normalize if not already done
if adata.raw is not None:
    print("    Using adata.raw for expression")
    expr_matrix = adata.raw
else:
    print("    Using adata.X directly")
    expr_matrix = adata

# Check key markers
for gene in ["TRPM5", "POU2F3", "KRT5", "KLF4", "CHAT"]:
    if gene in adata.var_names:
        # Get expression
        if adata.raw is not None:
            expr = np.array(adata.raw[:, gene].X.todense()).flatten() if hasattr(adata.raw[:, gene].X, 'todense') else np.array(adata.raw[:, gene].X).flatten()
        else:
            expr = np.array(adata[:, gene].X.todense()).flatten() if hasattr(adata[:, gene].X, 'todense') else np.array(adata[:, gene].X).flatten()
        nz = (expr > 0).sum()
        pct = nz / len(expr) * 100
        print(f"    {gene}: {nz} / {len(expr)} cells ({pct:.2f}%) express > 0")

# ============================================================
# 5. Map cell types using known annotations or marker scoring
# ============================================================
print(f"\n[5] Cell type annotation...")

# Find the cell type column
cell_type_col = None
for col in adata.obs.columns:
    if any(kw == col.lower() for kw in ['cell_type', 'celltype', 'cell.type', 'ct', 'cluster_name', 'annotation']):
        cell_type_col = col
        break
# Fallback: look for columns with 'cell' in name
if cell_type_col is None:
    for col in adata.obs.columns:
        if 'cell' in col.lower() and 'type' in col.lower():
            cell_type_col = col
            break

if cell_type_col:
    print(f"    Found cell type column: '{cell_type_col}'")
    ct_counts = adata.obs[cell_type_col].value_counts()
    print(f"    Cell types ({len(ct_counts)}):")
    for ct, count in ct_counts.items():
        pct = count / adata.n_obs * 100
        marker = " <- BRUSH/TUFT" if any(k in str(ct).lower() for k in ['tuft', 'brush', 'ionocyte']) else ""
        marker += " <- BASAL" if any(k in str(ct).lower() for k in ['basal']) else ""
        print(f"      {ct}: {count} cells ({pct:.1f}%){marker}")
else:
    print("    No cell type annotation found. Using marker-based scoring.")
    # Use scanpy to score cells for tuft and basal signatures
    if "TRPM5" in adata.var_names and "POU2F3" in adata.var_names:
        sc.tl.score_genes(adata, gene_list=["TRPM5", "POU2F3", "GNAT3"], score_name="tuft_score")
    if "KRT5" in adata.var_names and "TP63" in adata.var_names:
        sc.tl.score_genes(adata, gene_list=["KRT5", "TP63", "KRT14"], score_name="basal_score")
    if "KLF4" in adata.var_names and "KRT13" in adata.var_names:
        sc.tl.score_genes(adata, gene_list=["KLF4", "KRT13", "SOX2"], score_name="squamous_score")
    print("    Created scores: tuft_score, basal_score, squamous_score")

# ============================================================
# 6. Expression heatmap of key genes across cell types
# ============================================================
print(f"\n[6] Generating expression heatmap...")

if cell_type_col:
    # Subset to genes of interest that are available
    plot_genes = [g for g in ALL_GENES if g in adata.var_names]
    if len(plot_genes) > 3:
        # Get mean expression per cell type
        sc.tl.rank_genes_groups(adata, groupby=cell_type_col, n_genes=5, method='wilcoxon')

        # Dotplot of key genes across cell types
        fig, ax = plt.subplots(figsize=(16, 8))
        sc.pl.dotplot(adata, var_names=plot_genes, groupby=cell_type_col,
                      standard_scale='var', ax=ax, show=False)
        ax.set_title('Key Gene Expression Across Airway Cell Types\n'
                     '(GSE193816: Endobronchial Brush, Allergic Asthma Study)',
                     fontsize=12, fontweight='bold')
        plt.tight_layout()
        plt.savefig("D:/ai agent/first-cc/scrna_analysis/expression_dotplot.png",
                    dpi=150, bbox_inches='tight')
        print("    Saved: expression_dotplot.png")
else:
    print("    Skipping (no cell type column)")

# ============================================================
# 7. Violin plots for key genes
# ============================================================
print(f"\n[7] Generating violin plots...")
key_vln_genes = [g for g in ["TRPM5", "POU2F3", "CHAT", "KLF4", "KRT13", "KRT5", "CHRNA3"] if g in adata.var_names]
if cell_type_col and len(key_vln_genes) > 0:
    fig, axes = plt.subplots(1, len(key_vln_genes), figsize=(4*len(key_vln_genes), 6))
    if len(key_vln_genes) == 1:
        axes = [axes]
    for i, gene in enumerate(key_vln_genes):
        sc.pl.violin(adata, keys=gene, groupby=cell_type_col, ax=axes[i],
                     rotation=90, show=False)
        axes[i].set_title(gene, fontsize=10, fontweight='bold')
    plt.suptitle('Gene Expression Across Airway Cell Types\n(GSE193816)',
                 fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.savefig("D:/ai agent/first-cc/scrna_analysis/violin_plots.png",
                dpi=150, bbox_inches='tight')
    print("    Saved: violin_plots.png")

# ============================================================
# 8. LIANA Cell-Cell Communication Analysis
# ============================================================
print(f"\n[8] Running LIANA cell-cell communication analysis...")

try:
    import liana as li

    # Prepare data
    # LIANA needs: adata with cell type labels, raw counts, and condition info
    if cell_type_col:
        print(f"    Running LIANA on cell types from '{cell_type_col}'...")

        # Step 1: Prep - ensure raw counts are accessible
        if adata.raw is not None:
            li_adata = adata.raw.to_adata()
        else:
            li_adata = adata.copy()

        # Step 2: Run all methods
        # Check if we have condition/groups
        condition_col = None
        for col in adata.obs.columns:
            if any(kw in col.lower() for kw in ['condition', 'group', 'disease', 'status']):
                if adata.obs[col].nunique() >= 2:
                    condition_col = col
                    break

        if condition_col:
            print(f"    Condition column: '{condition_col}' ({sorted(adata.obs[condition_col].unique().astype(str))})")
            conditions = sorted(adata.obs[condition_col].unique().astype(str))
        else:
            conditions = ['all']
            print(f"    No condition column found, running on all cells")

        # For each condition, run LIANA
        all_lr_results = {}
        for cond in conditions:
            if condition_col and cond != 'all':
                cond_mask = adata.obs[condition_col].astype(str) == cond
                cond_adata = adata[cond_mask].copy()
                if cond_adata.n_obs < 100:
                    print(f"    Skipping '{cond}' (< 100 cells)")
                    continue
            else:
                cond_adata = adata.copy()

            print(f"\n    Condition: {cond} ({cond_adata.n_obs} cells)")

            # Run LIANA
            try:
                li.mt.rank_aggregate(
                    cond_adata,
                    groupby=cell_type_col,
                    resource_name='consensus',  # Uses multiple resources (CellPhoneDB, CellChat, etc.)
                    expr_prop=0.1,  # Min fraction of cells expressing the gene
                    n_perms=100,
                    verbose=True,
                    use_raw=False,
                )

                # Store results
                if hasattr(cond_adata, 'uns') and 'liana_res' in cond_adata.uns:
                    lr_res = cond_adata.uns['liana_res']
                    all_lr_results[cond] = lr_res
                    print(f"    Found {len(lr_res)} ligand-receptor interactions")

                    # Filter for tuft->basal specific interactions
                    tuft_types = [ct for ct in ct_counts.index if any(
                        k in str(ct).lower() for k in ['tuft', 'brush', 'ionocyte'])]
                    basal_types = [ct for ct in ct_counts.index if any(
                        k in str(ct).lower() for k in ['basal'])]

                    if tuft_types and basal_types:
                        tuft_to_basal = lr_res[
                            lr_res['source'].isin(tuft_types) &
                            lr_res['target'].isin(basal_types)
                        ]
                        print(f"    Tuft->Basal interactions: {len(tuft_to_basal)}")
                        if len(tuft_to_basal) > 0:
                            print(f"\n    Top tuft->basal LR pairs:")
                            top = tuft_to_basal.nsmallest(10, 'magnitude_rank')
                            for _, row in top.iterrows():
                                print(f"      {row['source']} -> {row['target']}: "
                                      f"{row['ligand_complex']} -> {row['receptor_complex']} "
                                      f"(rank={row['magnitude_rank']:.2f}, "
                                      f"p={row.get('pvalue', '?')})")

                            # Highlight cholinergic-specific pairs
                            cholinergic_pairs = tuft_to_basal[
                                tuft_to_basal['ligand_complex'].str.contains('ACH|CHAT|SLC18', case=False, na=False) |
                                tuft_to_basal['receptor_complex'].str.contains('CHRN|CHRM|nAChR|AChR', case=False, na=False)
                            ]
                            if len(cholinergic_pairs) > 0:
                                print(f"\n    *** CHOLINERGIC tuft->basal pairs ({len(cholinergic_pairs)}): ***")
                                for _, row in cholinergic_pairs.iterrows():
                                    print(f"      {row['ligand_complex']} -> {row['receptor_complex']} "
                                          f"(rank={row['magnitude_rank']:.2f})")

            except Exception as e:
                print(f"    LIANA error for '{cond}': {e}")
                import traceback
                traceback.print_exc()

    else:
        print("    Skipping LIANA: no cell type column found")

except ImportError as e:
    print(f"    LIANA import error: {e}")
    print("    Will generate manual LR dotplot instead.")

    # Manual LR analysis: compute mean expression of ligand/receptor pairs
    if cell_type_col:
        print("\n    Running manual LR analysis...")

        # Define LR pairs of interest
        lr_pairs = [
            # (ligand, receptor, pathway_name)
            ("CHAT", "CHRNA3", "ACh-nAChR α3"),
            ("CHAT", "CHRNA5", "ACh-nAChR α5"),
            ("CHAT", "CHRNA7", "ACh-nAChR α7"),
            ("CHAT", "CHRM3", "ACh-mAChR M3"),
            ("IL25", "IL17RB", "IL-25"),
            ("LTC4S", "CYSLTR1", "CysLT"),
            ("LTC4S", "CYSLTR2", "CysLT-2"),
        ]

        available_lr = []
        for lig, rec, name in lr_pairs:
            if lig in adata.var_names and rec in adata.var_names:
                available_lr.append((lig, rec, name))

        if available_lr:
            cell_types = sorted(adata.obs[cell_type_col].unique().astype(str))
            n_ct = len(cell_types)

            # Compute mean expression per cell type
            mean_expr = {}
            for gene in set([l for l, r, n in available_lr] + [r for l, r, n in available_lr]):
                if gene in adata.var_names:
                    gene_idx = list(adata.var_names).index(gene)
                    if hasattr(adata.X, 'todense'):
                        gene_expr = np.array(adata.X[:, gene_idx].todense()).flatten()
                    else:
                        gene_expr = np.array(adata.X[:, gene_idx]).flatten()
                    for ct in cell_types:
                        mask = adata.obs[cell_type_col].astype(str) == ct
                        mean_expr[(ct, gene)] = np.mean(gene_expr[mask]) if mask.sum() > 0 else 0

            # Build LR score matrix
            print(f"\n    Manual LR scores (ligand_expr * receptor_expr per cell type pair):")
            for lig, rec, name in available_lr:
                for src_ct in cell_types:
                    for tgt_ct in cell_types:
                        lig_exp = mean_expr.get((src_ct, lig), 0)
                        rec_exp = mean_expr.get((tgt_ct, rec), 0)
                        score = lig_exp * rec_exp
                        if score > 0.01:
                            print(f"      {name}: {src_ct} -> {tgt_ct}: {score:.4f}")

# ============================================================
# 9. UMAP visualization
# ============================================================
print(f"\n[9] Generating UMAP...")
if cell_type_col:
    fig, axes = plt.subplots(1, 2, figsize=(18, 7))

    # UMAP colored by cell type
    sc.pl.umap(adata, color=cell_type_col, ax=axes[0], show=False,
               title=f'Cell Types ({cell_type_col})', legend_loc='right margin')

    # UMAP of key genes
    umap_genes = [g for g in ["TRPM5", "KLF4", "KRT5"] if g in adata.var_names]
    if umap_genes:
        sc.pl.umap(adata, color=umap_genes[0], ax=axes[1], show=False,
                   title=f'{umap_genes[0]} expression', cmap='Reds')

    plt.tight_layout()
    plt.savefig("D:/ai agent/first-cc/scrna_analysis/umap_celltypes.png",
                dpi=150, bbox_inches='tight')
    print("    Saved: umap_celltypes.png")
else:
    # Try to find UMAP coordinates
    has_umap = all(k in adata.obsm for k in ['X_umap'])
    if has_umap:
        fig, ax = plt.subplots(figsize=(10, 8))
        sc.pl.umap(adata, ax=ax, show=False)
        plt.savefig("D:/ai agent/first-cc/scrna_analysis/umap.png",
                    dpi=150, bbox_inches='tight')
        print("    Saved: umap.png")

# ============================================================
# 10. Summary Report
# ============================================================
print("\n" + "=" * 70)
print("ANALYSIS SUMMARY")
print("=" * 70)

# Count tuft and basal cells
if cell_type_col:
    ct_counts = adata.obs[cell_type_col].value_counts()
    tuft_cts = [ct for ct in ct_counts.index if any(
        k in str(ct).lower() for k in ['tuft', 'brush'])]
    basal_cts = [ct for ct in ct_counts.index if any(
        k in str(ct).lower() for k in ['basal'])]

    tuft_n = sum(ct_counts[ct] for ct in tuft_cts)
    basal_n = sum(ct_counts[ct] for ct in basal_cts)

    print(f"\n  Tuft/Brush cells: {tuft_n} ({tuft_n/adata.n_obs*100:.2f}%)")
    for ct in tuft_cts:
        print(f"    - {ct}: {ct_counts[ct]} cells")

    print(f"\n  Basal cells: {basal_n} ({basal_n/adata.n_obs*100:.2f}%)")
    for ct in basal_cts:
        print(f"    - {ct}: {ct_counts[ct]} cells")

    # Check expression
    print("\n  Gene expression summary (mean in each population):")
    for gene in available_genes[:15]:
        if gene in adata.var_names:
            gene_idx = list(adata.var_names).index(gene)
            if hasattr(adata.X, 'todense'):
                gene_expr = np.array(adata.X[:, gene_idx].todense()).flatten()
            else:
                gene_expr = np.array(adata.X[:, gene_idx]).flatten()

            if tuft_cts:
                tuft_mask = adata.obs[cell_type_col].astype(str).isin(tuft_cts)
                tuft_mean = np.mean(gene_expr[tuft_mask]) if tuft_mask.sum() > 0 else 0
                tuft_pct = np.mean(gene_expr[tuft_mask] > 0) * 100 if tuft_mask.sum() > 0 else 0
            else:
                tuft_mean, tuft_pct = 0, 0

            if basal_cts:
                basal_mask = adata.obs[cell_type_col].astype(str).isin(basal_cts)
                basal_mean = np.mean(gene_expr[basal_mask]) if basal_mask.sum() > 0 else 0
                basal_pct = np.mean(gene_expr[basal_mask] > 0) * 100 if basal_mask.sum() > 0 else 0
            else:
                basal_mean, basal_pct = 0, 0

            tuft_str = f"Tuft={tuft_mean:.3f} ({tuft_pct:.0f}%)" if tuft_n > 0 else "N/A"
            basal_str = f"Basal={basal_mean:.3f} ({basal_pct:.0f}%)" if basal_n > 0 else "N/A"
            print(f"    {gene:12s}: {tuft_str:25s} | {basal_str}")

print(f"\nAnalysis complete.")
print("=" * 70)
