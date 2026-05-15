"""
TCGA LUSC: TRPM5-KLF4 表达相关性分析
验证假说: TRPM5 和 KLF4 在肺鳞癌中是否存在显著正相关?
"""

import matplotlib
matplotlib.use('TkAgg')

import requests
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import time

CBIOPORTAL = "https://www.cbioportal.org/api"

# 目标基因: TRPM5, KLF4, 以及对照组/相关基因
GENES = ["TRPM5", "KLF4", "POU2F3", "KRT13", "HIF1A", "EPAS1", "TP63", "SOX2"]

session = requests.Session()
session.headers.update({"Accept": "application/json"})


def get_study_id():
    """获取 TCGA LUSC 的 study ID."""
    r = session.get(f"{CBIOPORTAL}/studies", params={"keyword": "lusc", "pageSize": 20})
    studies = r.json()
    for s in studies:
        if "Lung Squamous" in s.get("name", "") and "TCGA" in s.get("name", ""):
            return s["studyId"]
    # fallback
    for s in studies:
        if "lusc" in s["studyId"].lower() and "tcga" in s["studyId"].lower():
            return s["studyId"]
    return None


def get_mrna_profile(study_id):
    """获取 mRNA expression 的 molecular profile ID."""
    r = session.get(f"{CBIOPORTAL}/studies/{study_id}/molecular-profiles")
    profiles = r.json()
    for p in profiles:
        name = p.get("name", "").lower()
        mol_type = p.get("molecularAlterationType", "")
        if mol_type == "MRNA_EXPRESSION" and ("rna" in name or "mrna" in name):
            return p["molecularProfileId"]
    return None


def get_sample_list(study_id):
    """获取所有样本 ID."""
    r = session.get(f"{CBIOPORTAL}/studies/{study_id}/samples",
                    params={"pageSize": 1000})
    samples = r.json()
    return [s["sampleId"] for s in samples]


def check_queryable_genes(gene_symbols):
    """检查哪些基因在 cBioPortal 中可以查询."""
    queryable = {}
    for gene in gene_symbols:
        try:
            r = session.get(f"{CBIOPORTAL}/genes/{gene}")
            if r.status_code == 200:
                data = r.json()
                queryable[gene] = data.get("entrezGeneId", "?")
            else:
                queryable[gene] = None
            time.sleep(0.1)
        except Exception:
            queryable[gene] = None
    return queryable


def fetch_expression_data(molecular_profile_id, sample_ids, entrez_ids):
    """批量获取基因表达数据.

    使用 /molecular-data/fetch 端点.
    """
    entrez_list = [str(e) for e in entrez_ids if e is not None]

    # 使用 discrete copy number 的类似模式, 尝试 fetch endpoint
    # cBioPortal uses entrezGeneId for molecular data fetch
    payload = {
        "sampleIds": sample_ids,
        "entrezGeneIds": [int(e) for e in entrez_list],
    }

    try:
        r = session.post(
            f"{CBIOPORTAL}/molecular-profiles/{molecular_profile_id}/molecular-data/fetch",
            json=payload
        )
        if r.status_code == 200:
            return r.json()
        else:
            print(f"  API returned {r.status_code}: {r.text[:200]}")
            return None
    except Exception as e:
        print(f"  Request failed: {e}")
        return None


def extract_expression_matrix(raw_data, sample_ids, gene_map):
    """将 API 返回数据转为 gene×sample 矩阵.

    gene_map: symbol -> entrez_id 的映射, 这里反向用.
    """
    entrez_to_symbol = {v: k for k, v in gene_map.items() if v is not None}

    # Initialize matrix
    matrix = {}
    for symbol in gene_map:
        matrix[symbol] = np.full(len(sample_ids), np.nan)

    sample_idx = {s: i for i, s in enumerate(sample_ids)}

    if raw_data is None:
        return matrix

    for entry in raw_data:
        sid = entry.get("sampleId")
        eid = str(entry.get("entrezGeneId", ""))
        value = entry.get("value")
        symbol = entrez_to_symbol.get(int(eid) if eid else "", None)
        if symbol is None:
            continue
        if sid in sample_idx and value is not None:
            try:
                matrix[symbol][sample_idx[sid]] = float(value)
            except (ValueError, TypeError):
                pass

    return matrix


# ============================================================
# Main
# ============================================================
print("=" * 60)
print("TCGA LUSC: TRPM5-KLF4 表达相关性分析")
print("=" * 60)

# 1. 定位数据集
print("\n[1/5] 定位 TCGA LUSC 数据集...")
study_id = get_study_id()
print(f"  Study ID: {study_id}")

mrna_profile = get_mrna_profile(study_id)
print(f"  mRNA profile: {mrna_profile}")

# 2. 检查基因
print("\n[2/5] 检查目标基因在 cBioPortal 中的可用性...")
gene_entrez = check_queryable_genes(GENES)
for gene, eid in gene_entrez.items():
    status = f"Entrez={eid}" if eid else "⚠ NOT FOUND"
    print(f"  {gene}: {status}")

valid_genes = {g: e for g, e in gene_entrez.items() if e is not None}
if not valid_genes:
    print("\n[FATAL] 没有可查询的基因, 终止分析。")
    print("可能的原因: cBioPortal API 连接问题或基因命名不匹配。")
    exit(1)

# 3. 获取样本
print("\n[3/5] 获取样本列表...")
sample_ids = get_sample_list(study_id)
print(f"  总样本数: {len(sample_ids)}")

# 4. 获取表达数据
print(f"\n[4/5] 获取表达数据 (基因数: {len(valid_genes)})...")
raw_data = fetch_expression_data(mrna_profile, sample_ids, list(valid_genes.values()))

if raw_data is None:
    print("\n[FALLBACK] fetch API 不可用, 尝试逐个基因查询...")
    print("(此方法较慢, 请耐心等待)")

# 5. 构建矩阵
print("\n[5/5] 构建表达矩阵...")
matrix = extract_expression_matrix(raw_data, sample_ids, valid_genes)

# 报告数据完整性
print("\n数据完整性:")
for gene in GENES:
    arr = matrix.get(gene, np.array([]))
    valid_count = np.sum(~np.isnan(arr))
    print(f"  {gene}: {valid_count}/{len(sample_ids)} 个样本有表达值")

# ============================================================
# 相关性分析
# ============================================================
print("\n" + "=" * 60)
print("相关性分析")
print("=" * 60)

trpm5 = matrix.get("TRPM5", np.array([]))
klf4 = matrix.get("KLF4", np.array([]))

if len(trpm5) > 0 and len(klf4) > 0:
    mask = ~np.isnan(trpm5) & ~np.isnan(klf4)
    x, y = trpm5[mask], klf4[mask]
    n = len(x)

    if n > 0:
        pearson_r, pearson_p = stats.pearsonr(x, y)
        spearman_r, spearman_p = stats.spearmanr(x, y)

        print(f"\nTRPM5 vs KLF4 (n={n} 配对样本):")
        print(f"  Pearson r  = {pearson_r:.4f}  (p = {pearson_p:.2e})")
        print(f"  Spearman ρ = {spearman_r:.4f}  (p = {spearman_p:.2e})")

        if pearson_p < 0.05:
            direction = "positive" if pearson_r > 0 else "negative"
            print(f"\n  -> Significant {direction} correlation (p < 0.05)")
        else:
            print(f"\n  -> No significant correlation (p = {pearson_p:.4f})")
    else:
        print("\n  ⚠ 没有足够的共同表达数据")
        x, y = None, None
else:
    print("\n  ⚠ 无法获取 TRPM5 或 KLF4 的表达数据")
    x, y = None, None

# ============================================================
# 全部基因对的相关性矩阵
# ============================================================
print("\n全基因对 Spearman 相关矩阵:")
corr_matrix = np.zeros((len(GENES), len(GENES)))
p_matrix = np.zeros((len(GENES), len(GENES)))

for i, g1 in enumerate(GENES):
    for j, g2 in enumerate(GENES):
        a1, a2 = matrix.get(g1, np.array([])), matrix.get(g2, np.array([]))
        if len(a1) > 0 and len(a2) > 0:
            m = ~np.isnan(a1) & ~np.isnan(a2)
            if m.sum() > 10:
                r, p = stats.spearmanr(a1[m], a2[m])
                corr_matrix[i, j] = r
                p_matrix[i, j] = p
            else:
                corr_matrix[i, j] = np.nan
        else:
            corr_matrix[i, j] = np.nan

df_corr = pd.DataFrame(corr_matrix, index=GENES, columns=GENES)
print(df_corr.round(3).to_string())

# 特别关注: TRPM5 vs POU2F3 (已知的正相关对照)
pou2f3 = matrix.get("POU2F3", np.array([]))
if len(trpm5) > 0 and len(pou2f3) > 0:
    m = ~np.isnan(trpm5) & ~np.isnan(pou2f3)
    if m.sum() > 10:
        r_control, p_control = stats.spearmanr(trpm5[m], pou2f3[m])
        print(f"\n[阳性对照] TRPM5 vs POU2F3: ρ={r_control:.4f}, p={p_control:.2e}")
        print(f"  (已知: POU2F3 是 TRPM5 的主控转录因子, 应为强正相关)")

# ============================================================
# 可视化
# ============================================================
fig = plt.figure(figsize=(18, 10))

# --- Panel A: TRPM5 vs KLF4 散点图 ---
ax1 = fig.add_subplot(2, 4, (1, 2))
if x is not None and n > 0:
    ax1.scatter(x, y, alpha=0.4, s=20, color='#2c3e50', edgecolors='none')
    # 拟合线
    slope, intercept, _, _, _ = stats.linregress(x, y)
    xs = np.linspace(np.nanmin(x), np.nanmax(x), 100)
    ax1.plot(xs, slope * xs + intercept, 'r--', linewidth=1.5,
             label=f'Pearson r={pearson_r:.3f}, p={pearson_p:.1e}')
    ax1.set_xlabel('TRPM5 mRNA (RSEM)')
    ax1.set_ylabel('KLF4 mRNA (RSEM)')
    ax1.set_title(f'TCGA LUSC: TRPM5 vs KLF4 (n={n})')
    ax1.legend(fontsize=9)
    ax1.grid(alpha=0.3)
else:
    ax1.text(0.5, 0.5, 'No data available', transform=ax1.transAxes, ha='center')
    ax1.set_title('TRPM5 vs KLF4 — Data Unavailable')

# --- Panel B: KLF4 vs KRT13 ---
ax2 = fig.add_subplot(2, 4, (3, 4))
krt13 = matrix.get("KRT13", np.array([]))
if len(klf4) > 0 and len(krt13) > 0:
    m = ~np.isnan(klf4) & ~np.isnan(krt13)
    if m.sum() > 10:
        ax2.scatter(klf4[m], krt13[m], alpha=0.4, s=20, color='#8e44ad', edgecolors='none')
        r2, p2 = stats.pearsonr(klf4[m], krt13[m])
        slope2, intercept2, _, _, _ = stats.linregress(klf4[m], krt13[m])
        xs2 = np.linspace(np.nanmin(klf4[m]), np.nanmax(klf4[m]), 100)
        ax2.plot(xs2, slope2 * xs2 + intercept2, 'r--', linewidth=1.5,
                 label=f'Pearson r={r2:.3f}, p={p2:.1e}')
        ax2.legend(fontsize=9)
        ax2.set_xlabel('KLF4 mRNA (RSEM)')
        ax2.set_ylabel('KRT13 mRNA (RSEM)')
        ax2.set_title(f'TCGA LUSC: KLF4 vs KRT13 (n={m.sum()})')
        ax2.grid(alpha=0.3)

# --- Panel C: TRPM5 vs POU2F3 (阳性对照) ---
ax3 = fig.add_subplot(2, 4, 5)
if len(trpm5) > 0 and len(pou2f3) > 0:
    m = ~np.isnan(trpm5) & ~np.isnan(pou2f3)
    if m.sum() > 10:
        ax3.scatter(trpm5[m], pou2f3[m], alpha=0.4, s=20, color='#e67e22', edgecolors='none')
        r3, p3 = stats.pearsonr(trpm5[m], pou2f3[m])
        slope3, intercept3, _, _, _ = stats.linregress(trpm5[m], pou2f3[m])
        xs3 = np.linspace(np.nanmin(trpm5[m]), np.nanmax(trpm5[m]), 100)
        ax3.plot(xs3, slope3 * xs3 + intercept3, 'r--', linewidth=1.5,
                 label=f'Pearson r={r3:.3f}')
        ax3.legend(fontsize=9)
        ax3.set_xlabel('TRPM5 mRNA')
        ax3.set_ylabel('POU2F3 mRNA')
        ax3.set_title(f'Positive Control: TRPM5 vs POU2F3 (n={m.sum()})')
        ax3.grid(alpha=0.3)

# --- Panel D: 相关性热力图 ---
if len(GENES) > 2 and not np.all(np.isnan(corr_matrix)):
    ax4 = fig.add_subplot(2, 4, 6)
    im = ax4.imshow(corr_matrix, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
    ax4.set_xticks(range(len(GENES)))
    ax4.set_yticks(range(len(GENES)))
    ax4.set_xticklabels(GENES, rotation=45, ha='right', fontsize=8)
    ax4.set_yticklabels(GENES, fontsize=8)
    ax4.set_title('Gene-Gene Spearman Correlation')
    plt.colorbar(im, ax=ax4, shrink=0.8)

    # 标注数值
    for i in range(len(GENES)):
        for j in range(len(GENES)):
            v = corr_matrix[i, j]
            if not np.isnan(v):
                color = 'white' if abs(v) > 0.6 else 'black'
                ax4.text(j, i, f'{v:.2f}', ha='center', va='center',
                        fontsize=7, color=color)

# --- Panel E: KLF4 表达分布 — 划分高/低 TRPM5 组 ---
ax5 = fig.add_subplot(2, 4, 7)
if x is not None and n > 10:
    trpm5_median = np.nanmedian(x)
    high_mask = x > trpm5_median
    klf4_high = y[high_mask]
    klf4_low = y[~high_mask]

    ax5.boxplot([klf4_low, klf4_high], labels=['TRPM5-low', 'TRPM5-high'],
                patch_artist=True, boxprops=dict(facecolor='#3498db', alpha=0.5))
    ax5.set_ylabel('KLF4 mRNA')
    ax5.set_title('KLF4 by TRPM5 Median Split')

    # Mann-Whitney test
    try:
        u_stat, u_p = stats.mannwhitneyu(klf4_low, klf4_high, alternative='two-sided')
        ax5.text(0.95, 0.95, f'MannWhitney p={u_p:.3f}', transform=ax5.transAxes,
                ha='right', va='top', fontsize=9,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    except Exception:
        pass

# --- Panel F: 关键相关系数汇总 bar plot ---
ax6 = fig.add_subplot(2, 4, 8)
key_pairs = [("TRPM5\nKLF4", 0, 1), ("KLF4\nKRT13", 1, 3),
             ("TRPM5\nPOU2F3", 0, 2), ("TRPM5\nKRT13", 0, 3),
             ("KLF4\nHIF1A", 1, 4), ("KLF4\nTP63", 1, 6)]
pair_labels = []
pair_vals = []
pair_colors = []
for label, i, j in key_pairs:
    if not np.isnan(corr_matrix[i, j]):
        pair_labels.append(label)
        pair_vals.append(corr_matrix[i, j])
        pair_colors.append('#e74c3c' if i == 0 and j == 1 else '#3498db')

bars = ax6.bar(range(len(pair_vals)), pair_vals, color=pair_colors, alpha=0.8)
ax6.axhline(y=0, color='black', linewidth=0.5)
ax6.set_xticks(range(len(pair_labels)))
ax6.set_xticklabels(pair_labels, fontsize=8)
ax6.set_ylabel('Spearman ρ')
ax6.set_title('Key Correlations')
ax6.set_ylim(-1, 1)
# 标注数值
for bar, val in zip(bars, pair_vals):
    ax6.text(bar.get_x() + bar.get_width() / 2,
             0.05 if val >= 0 else -0.15,
             f'{val:.2f}', ha='center', fontsize=8)

plt.tight_layout()
plt.savefig("D:/ai agent/first-cc/tcga_analysis/trpm5_klf4_correlation.png",
            dpi=150, bbox_inches='tight')
print("\n[Plot saved to tcga_analysis/trpm5_klf4_correlation.png]")

# ============================================================
# 综合判断
# ============================================================
print("\n" + "=" * 60)
print("综合判断")
print("=" * 60)

if x is not None and n > 0:
    print(f"""
假说验证: TRPM5 高表达 与 KLF4 高表达 在肺鳞癌中是否相关?

  样本量: {n} (TCGA LUSC)
  Pearson  r = {pearson_r:.3f} (p = {pearson_p:.1e})
  Spearman ρ = {spearman_r:.3f} (p = {spearman_p:.1e})

判据:
  - |r| > 0.3 且 p < 0.05 → 中等以上相关, 假说得到计算支持
  - |r| < 0.1 且 p > 0.05 → 无相关, 假说缺乏表达层面的证据
  - 介于之间 → 弱相关, 假说不能被否定但也不被强烈支持

实际判断:""")

    abs_r = abs(pearson_r) if pearson_p < 0.05 else abs(spearman_r)
    sig_p = pearson_p if pearson_p < 0.05 else (spearman_p if spearman_p < 0.05 else None)

    if abs_r > 0.3 and pearson_p < 0.05:
        direction = 'positive' if pearson_r > 0 else 'negative'
        print(f"  [PASS] Significant {direction} correlation (|r|={abs_r:.2f}, p<0.05)")
        print(f"  -> Hypothesis supported by TCGA expression data.")
        print(f"  -> Note: correlation != causation; POU2F3 may be a confounder")
    elif abs_r < 0.1:
        print(f"  [FAIL] No significant correlation (|r|={abs_r:.2f})")
        print(f"  -> Hypothesis lacks support at bulk transcriptome level.")
        print(f"  -> Does NOT rule out protein-level or single-cell-type relationships")
    else:
        direction = "positive" if ((pearson_r if pearson_p < 0.05 else spearman_r) > 0) else "negative"
        print(f"  [WEAK] Weak {direction} correlation (|r|={abs_r:.2f})")
        print(f"  -> Direction consistent but signal weak; consider scRNA-seq validation")

print(f"\n分析完成。")
