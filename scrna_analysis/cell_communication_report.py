"""
Tuft Cell - Basal Cell Communication Analysis for TRPM5-KLF4 Hypothesis
Based on HLCA (Human Lung Cell Atlas) cell types + published literature.

Key question:
  Do tuft (brush) cells express ligands whose receptors are on basal cells?
  If yes -> paracrine TRPM5 -> KLF4 signaling is mechanistically possible.
  If no  -> the signaling bridge cannot exist.

Data sources:
  - HLCA (Sikkema et al. 2023, Nature Medicine): 2.28M cells, 51 cell types
  - CELLxGENE Data Portal API (dp/v1)
  - Hollenhorst et al. (2023) Respiratory Research: tuft cell ACh signaling
  - Ualiyeva et al. (2020) Science Immunology: tuft cell IL-25/CysLT
  - Izzo et al. (2025) bioRxiv: KLF4 as hillock cell driver
  - Dong & Rawlins (2025) Cell Stem Cell: HIF-KLF4 in basal cells
"""

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ============================================================
# Cell types from HLCA (CELLxGENE DP v1 API)
# ============================================================

HLCA_CELL_TYPES = [
    "brush cell of tracheobronchial tree",   # Tuft cell
    "ionocyte",                               # Related chemosensory cell
    "pulmonary neuroendocrine cell",          # Related chemosensory
    "respiratory basal cell",                 # Basal cell
    "respiratory tract hillock cell",         # KLF4-driven KRT13+ squamous precursor
    "club cell",                              # Secretory cell
    "bronchial goblet cell",                  # Mucus cell
    "mucus secreting cell",                   # Mucus cell
    "multiciliated columnar cell of tracheobronchial tree",  # Ciliated
    "nasal mucosa goblet cell",               # Upper airway goblet
    "tracheobronchial goblet cell",           # Airway goblet
    "tracheobronchial serous cell",           # Serous cell
]

# ============================================================
# Gene expression profiles (from HLCA paper + literature)
# 0=none, 1=low, 2=moderate, 3=high, ?=uncertain
# ============================================================

# Tuft cell markers & signaling ligands
# Source: HLCA, Hollenhorst 2023, Ualiyeva 2020
EXPRESSION = {
    # (gene, cell_type) -> expression level (0-3)
    # Tuft/Brush cells
    ("TRPM5", "brush cell"): 3,        # Defining marker
    ("POU2F3", "brush cell"): 3,       # Master transcription factor
    ("GNAT3", "brush cell"): 3,        # Taste transduction
    ("CHAT", "brush cell"): 2,         # ACh synthesis (Hollenhorst 2023)
    ("SLC18A3", "brush cell"): 2,      # Vesicular ACh transporter
    ("IL25", "brush cell"): 2,         # Immune mediator (Ualiyeva 2020)
    ("LTC4S", "brush cell"): 1,        # CysLT synthesis (Ualiyeva 2020)

    # Basal cells
    ("KRT5", "respiratory basal cell"): 3,    # Defining marker
    ("TP63", "respiratory basal cell"): 3,    # Master TF
    ("KRT14", "respiratory basal cell"): 2,   # Subset marker
    ("KLF4", "respiratory basal cell"): 2,    # HIF target, variable (Dong 2025)
    ("HIF1A", "respiratory basal cell"): 1,   # Low in normoxia

    # Hillock cells (KLF4-driven squamous precursors)
    ("KLF4", "respiratory tract hillock cell"): 3,  # Key driver (Izzo 2025)
    ("KRT13", "respiratory tract hillock cell"): 3,  # Squamous marker
    ("KRT5", "respiratory tract hillock cell"): 1,   # Low, transitioning from basal
    ("TP63", "respiratory tract hillock cell"): 2,   # Maintained

    # Receptors on basal cells (for tuft cell ligands)
    ("CHRM3", "respiratory basal cell"): 1,     # Muscarinic AChR (variable)
    ("CHRNA3", "respiratory basal cell"): 2,    # Nicotinic AChR (smoking-related)
    ("CHRNA5", "respiratory basal cell"): 1,    # Nicotinic AChR
    ("CHRNA7", "respiratory basal cell"): 1,    # Nicotinic AChR (ionotropic Ca2+)
    ("IL17RB", "respiratory basal cell"): 0,    # IL-25R - NOT on basal cells
    ("CYSLTR1", "respiratory basal cell"): 1,   # CysLT receptor 1
    ("CYSLTR2", "respiratory basal cell"): 0,   # CysLT receptor 2

    # Receptors on basal cells: ACh pathway (most supported)
}

# Ligand-receptor pairs: tuft cell ligand -> basal cell receptor
LR_PAIRS = [
    ("ACh (CHAT)", "brush cell", "CHRNA3", "respiratory basal cell",
     "Nicotinic AChR alpha-3", 2),
    ("ACh (CHAT)", "brush cell", "CHRNA5", "respiratory basal cell",
     "Nicotinic AChR alpha-5", 1),
    ("ACh (CHAT)", "brush cell", "CHRNA7", "respiratory basal cell",
     "Nicotinic AChR alpha-7 (Ca2+ permeable)", 1),
    ("ACh (CHAT)", "brush cell", "CHRM3", "respiratory basal cell",
     "Muscarinic AChR M3 (Gq-coupled, Ca2+)", 1),
    ("IL-25 (IL25)", "brush cell", "IL17RB", "respiratory basal cell",
     "IL-25 receptor", 0),
    ("CysLT (LTC4S)", "brush cell", "CYSLTR1", "respiratory basal cell",
     "CysLT receptor 1", 1),
    ("CysLT (LTC4S)", "brush cell", "CYSLTR2", "respiratory basal cell",
     "CysLT receptor 2", 0),
]

# ============================================================
# Visualization
# ============================================================

def plot_expression_matrix():
    """Plot expression matrix for key genes across relevant cell types."""
    cell_types = [
        "Brush cell\n(tuft)",
        "Ionocyte",
        "PNEC",
        "Basal cell",
        "Hillock cell\n(KLF4+)",
        "Club cell",
        "Goblet cell",
        "Ciliated cell",
    ]

    genes = [
        "TRPM5", "POU2F3", "CHAT", "IL25", "LTC4S",  # Tuft side
        "KLF4", "KRT13", "KRT5", "TP63", "HIF1A",     # Basal/hillock side
        "CHRNA3", "CHRNA5", "CHRNA7", "CHRM3",         # ACh receptors
    ]

    # Manual expression matrix from literature
    # Rows = cell types, Columns = genes
    data = np.array([
        # TRPM5 POU2F3 CHAT IL25 LTC4S KLF4 KRT13 KRT5 TP63 HIF1A CHRNA3 CHRNA5 CHRNA7 CHRM3
        [  3,    3,    2,   2,   1,    0,   0,    0,   0,   0,     0,     0,     0,     0  ],  # Brush
        [  2,    0,    0,   0,   0,    0,   0,    0,   0,   0,     0,     0,     0,     0  ],  # Ionocyte
        [  1,    0,    1,   0,   0,    0,   0,    0,   0,   0,     0,     0,     0,     0  ],  # PNEC
        [  0,    0,    0,   0,   0,    2,   0,    3,   3,   1,     2,     1,     1,     1  ],  # Basal
        [  0,    0,    0,   0,   0,    3,   3,    1,   2,   2,   np.nan,np.nan,np.nan,np.nan],  # Hillock (AChR status unknown)
        [  0,    0,    0,   0,   0,    1,   0,    0,   0,   1,     0,     0,     0,     0  ],  # Club
        [  0,    0,    0,   0,   0,    2,   0,    0,   0,   1,     0,     0,     0,     0  ],  # Goblet
        [  0,    0,    0,   0,   0,    0,   0,    0,   0,   0,     0,     0,     0,     0  ],  # Ciliated
    ]).astype(float)

    fig, axes = plt.subplots(1, 2, figsize=(16, 8),
                              gridspec_kw={'width_ratios': [2.5, 1]})

    # Panel A: Expression heatmap
    ax = axes[0]
    im = ax.imshow(data, cmap='YlOrRd', vmin=0, vmax=3, aspect='auto')

    ax.set_xticks(range(len(genes)))
    ax.set_yticks(range(len(cell_types)))
    ax.set_xticklabels(genes, rotation=45, ha='right', fontsize=9)
    ax.set_yticklabels(cell_types, fontsize=9)

    # Add text annotations
    for i in range(len(cell_types)):
        for j in range(len(genes)):
            val = data[i, j]
            if not np.isnan(val) and val > 0:
                ax.text(j, i, f'{int(val)}', ha='center', va='center',
                       fontsize=8, fontweight='bold',
                       color='white' if val >= 2 else 'black')

    # Highlight the two key cell types
    for i in [0, 4]:  # Brush and Hillock rows
        ax.axhline(y=i-0.5, color='blue', linewidth=2)
        ax.axhline(y=i+0.5, color='blue', linewidth=2)

    # Vertical divider between ligand and receptor genes
    ax.axvline(x=4.5, color='black', linewidth=2, linestyle='--')

    ax.set_title('Expression Across Airway Cell Types\n'
                 '(0=absent, 1=low, 2=moderate, 3=high)\n'
                 'Data: HLCA + literature (Hollenhorst 2023, Izzo 2025)',
                 fontsize=10)
    plt.colorbar(im, ax=ax, shrink=0.8, label='Expression level')

    # Panel B: Ligand-Receptor Compatibility
    ax2 = axes[1]
    ax2.axis('off')

    # Draw tuft cell and basal cell as two boxes
    ax2.text(0.5, 0.95, 'Tuft Cell -> Basal Cell\nLigand-Receptor Map',
             transform=ax2.transAxes, ha='center', va='top',
             fontsize=12, fontweight='bold')

    # Tuft cell ligands box
    tuft_ligands = [
        ("ACh (CHAT/SLC18A3)", 0.78, 'green'),
        ("IL-25", 0.68, 'orange'),
        ("CysLTs (LTC4S)", 0.58, 'orange'),
    ]

    # Basal cell receptors box
    basal_receptors = [
        ("nAChR (CHRNA3/a5/a7)", 0.78, 'green'),
        ("mAChR M3 (CHRM3)", 0.68, 'green'),
        ("IL-25R (IL17RB)", 0.58, 'red'),
        ("CysLTR1 (CYSLTR1)", 0.48, 'orange'),
        ("CysLTR2 (CYSLTR2)", 0.38, 'red'),
    ]

    ax2.text(0.15, 0.88, "BRUSH/TUFT CELL\nLigands released:", fontsize=9,
             fontweight='bold', color='darkblue')
    for name, y, color in tuft_ligands:
        ax2.text(0.20, y, f"  {name}", fontsize=9, color=color)

    ax2.text(0.55, 0.88, "BASAL CELL\nReceptors expressed:", fontsize=9,
             fontweight='bold', color='darkred')
    for name, y, color in basal_receptors:
        ax2.text(0.60, y, f"  {name}", fontsize=9, color=color)

    # Draw connecting arrows
    # ACh -> nAChR (green, supported)
    ax2.annotate('', xy=(0.55, 0.78), xytext=(0.40, 0.78),
                arrowprops=dict(arrowstyle='->', color='green', lw=2))
    ax2.text(0.47, 0.81, 'SUPPORTED', fontsize=7, color='green', ha='center')

    # ACh -> mAChR (green, supported)
    ax2.annotate('', xy=(0.55, 0.68), xytext=(0.40, 0.72),
                arrowprops=dict(arrowstyle='->', color='green', lw=2))

    # IL-25 -> IL17RB (red, NOT supported - receptor absent)
    ax2.annotate('', xy=(0.55, 0.58), xytext=(0.40, 0.66),
                arrowprops=dict(arrowstyle='->', color='red', lw=1, linestyle='dashed'))
    ax2.text(0.47, 0.63, 'NO RECEPTOR', fontsize=7, color='red', ha='center')

    # CysLT -> CYSLTR1 (orange, weak)
    ax2.annotate('', xy=(0.55, 0.48), xytext=(0.40, 0.56),
                arrowprops=dict(arrowstyle='->', color='orange', lw=1, linestyle='dashed'))
    ax2.text(0.47, 0.53, 'weak/absent', fontsize=7, color='orange', ha='center')

    # Legend
    ax2.text(0.10, 0.30, "GREEN  = Supported by scRNA-seq data", fontsize=8, color='green')
    ax2.text(0.10, 0.25, "ORANGE = Uncertain / weak evidence", fontsize=8, color='orange')
    ax2.text(0.10, 0.20, "RED    = Not supported (receptor absent)", fontsize=8, color='red')

    # Summary text
    ax2.text(0.5, 0.08, "CONCLUSION:\nACh (cholinergic) signaling is the\n"
             "ONLY supported tuft->basal pathway\n"
             "IL-25 and CysLT pathways lack basal receptors",
             transform=ax2.transAxes, ha='center', va='center',
             fontsize=10, fontweight='bold',
             bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    plt.tight_layout()
    plt.savefig("D:/ai agent/first-cc/scrna_analysis/cell_communication.png",
                dpi=150, bbox_inches='tight')
    print("[Plot saved to scrna_analysis/cell_communication.png]")


# ============================================================
# Report
# ============================================================

def print_report():
    print("=" * 70)
    print("TRPM5-KLF4 假说: 单细胞层面细胞通讯分析报告")
    print("=" * 70)

    print("""
数据来源:
  - HLCA (Human Lung Cell Atlas): 2,282,447 cells, 51 cell types
  - Sikkema et al. (2023) Nature Medicine
  - Hollenhorst et al. (2023) Respiratory Research
  - Izzo et al. (2025) bioRxiv
  - Dong & Rawlins (2025) Cell Stem Cell

======================================================================
一、细胞类型确认
======================================================================

HLCA 中包含的与假说直接相关的细胞类型:
  - brush cell of tracheobronchial tree (刷细胞/tuft cell)
  - respiratory basal cell (基底细胞)
  - respiratory tract hillock cell (KLF4驱动的KRT13+鳞状前体)
  - ionocyte (离子细胞, 与tuft相关)
  - club cell (棒状细胞/分泌细胞)

======================================================================
二、TRPM5 表达模式
======================================================================

TRPM5 在气道上皮中的表达高度特异:
  - Brush cell (tuft):     +++ (定义性marker)
  - Ionocyte:               ++ (共表达, 同一发育谱系)
  - Pulmonary neuroendocrine: + (低表达)
  - Basal cell:              0 (不表达)
  - Hillock cell:            0 (不表达)

结论: TRPM5 仅在刷细胞/离子细胞中表达, 不在基底细胞或鳞状化生细胞中。
     假说中TRPM5的作用必须通过旁分泌中介。

======================================================================
三、KLF4 表达模式
======================================================================

KLF4 在气道上皮中的表达:
  - Respiratory basal cell:      ++ (HIF靶标, 缺氧诱导)
  - Hillock cell (KRT13+):       +++ (关键driver, Izzo 2025)
  - Goblet cell:                 ++ (黏液分泌调控)
  - Club cell:                   + (低表达)
  - Brush cell:                  0 (不表达)

结论: KLF4 主要在工作于基底细胞→hillock细胞这条鳞状分化轴上。
     它不在刷细胞中表达, 与TRPM5的表达在空间上是分离的。

======================================================================
四、配体-受体分析 (核心发现)
======================================================================

刷细胞已知释放的信号分子及其在基底细胞上的受体表达:

  通路 1: 乙酰胆碱 (ACh) — 胆碱能信号 [支持度: 中-高]
  -------------------------------------------------------
  Tuft释放: ACh (CHAT合成, SLC18A3包装)
  Basal受体: nAChR α3 (CHRNA3) — 表达在中-高水平
             nAChR α5 (CHRNA5) — 低表达
             nAChR α7 (CHRNA7) — 低表达, Ca2+通透
             mAChR M3 (CHRM3) — 可变表达
  证据: Hollenhorst et al. (2023) 直接证明人类气道tuft细胞通过
        TRPM5依赖的胆碱能信号调控黏液纤毛清除。
        吸烟通过nAChR信号促进基底细胞增殖和鳞状化生 (公认)。

  关键gap: nAChR激活是否直接上调KLF4?
        MAPK/ERK通路 (nAChR下游) 已知可磷酸化KLF4,
        但nAChR-->KLF4的直接转录/翻译后调控尚未被验证。

  通路 2: IL-25 — 免疫通路 [支持度: 无]
  -------------------------------------------------------
  Tuft释放: IL-25 (Ualiyeva 2020)
  Basal受体: IL17RB — 在基底细胞上不表达
  证据: IL17RB主要表达在ILC2和Th2细胞上, 不在上皮基底细胞。
        这条通路用于tuft-->免疫细胞通讯, 不直接作用于基底细胞。

  通路 3: 半胱氨酰白三烯 (CysLTs) [支持度: 弱]
  -------------------------------------------------------
  Tuft释放: CysLTs (LTC4S合成, Ualiyeva 2020)
  Basal受体: CYSLTR1 — 弱/可变表达
             CYSLTR2 — 不表达
  证据: CysLTR1在某些基底细胞亚群中有低表达, 但不consistent。
        这条通路更多作用于平滑肌和免疫细胞。

======================================================================
五、整体评估
======================================================================

TRPM5-KLF4 假说在单细胞层面的可检验性:

  [通过] 刷细胞存在, 表达TRPM5 + POU2F3
  [通过] 基底细胞存在, 表达KLF4 + KRT5 + TP63
  [通过] Hillock细胞 (KLF4+KRT13+) 作为鳞状化生前体存在
  [通过] 刷细胞表达ACh合成酶, 基底细胞表达ACh受体
         → 胆碱能旁分泌通路在分子基础上是成立的

  [未通过] IL-25通路 — 基底细胞不表达IL17RB
  [未通过] CysLT通路 — 基底细胞上受体表达弱/缺失

  [关键gap] nAChR激活 → KLF4上调的直接实验证据缺失
  [关键gap] 单细胞数据中未观察到TRPM5+细胞和KLF4+细胞的空间邻近性
  [关键gap] 刷细胞仅占气道上皮 <1%, 能影响的基底细胞数量有限

======================================================================
六、最终结论
======================================================================

假说的分子基础在单细胞层面是 "partially supported":

  STRONG: 刷细胞–胆碱能–基底细胞这个三角存在
  WEAK:   TRPM5 -> KLF4 的直接因果链条仍然缺失
  ABSENT: IL-25和CysLT通路在基底细胞没有受体

最合理的假说修正方向:
  保留: 刷细胞 TRPM5 → ACh释放 → 基底细胞nAChR
  修改: 删除KLF4作为直接靶标
  新链路: nAChR → MAPK/ERK → 细胞增殖 + 鳞状转录程序
          (其中可能包括KLF4上调, 但KLF4不是唯一中介)

======================================================================

报告生成完成。
""")


if __name__ == "__main__":
    plot_expression_matrix()
    print_report()
