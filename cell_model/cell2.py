"""
TRPM5-KLF4-KRT13 肺鳞癌化生动力学模型 (v2.0)
扩展内容:
  1. TRPM5 脱敏动力学 (报告建议 3.1.3)
  2. V_m → Ca²⁺ 中间信号层 (报告建议 3.1.1)
  3. HIF 正反馈恶性循环 (报告机制 4)
  4. 参数敏感性分析 (报告建议 3.1.2)
  5. 基于真实实验数据的参数校准 (报告第二节)
"""

import matplotlib
matplotlib.use('TkAgg')

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from model_params import MODEL_PARAMS as P, PARAM_RANGES


# ============================================================
# 模型 A: 基础模型 (原始结构, 参数校准自报告数据)
# ============================================================
class BaseModel:
    """与原始 cell2.py 结构一致, 但参数从 model_params 统一加载。"""

    def __init__(self):
        self.k_deg_klf4 = P["k_deg_klf4"]
        self.V_basal_klf4 = P["V_basal_klf4"]
        self.V_signal_max = P["V_signal_max"]
        self.n_signal = P["n_signal"]
        self.Kd_signal = P["Kd_signal"]
        self.V_auto = P["V_auto"]
        self.m_klf4 = P["m_klf4"]
        self.Kd_klf4_auto = P["Kd_klf4_auto"]
        self.k_krt13 = P["k_krt13"]
        self.k_deg_krt13 = P["k_deg_krt13"]

    def __call__(self, t, y, TRPM5_activity, stress_factor):
        KLF4, KRT13 = y
        effective = TRPM5_activity * stress_factor

        drive = self.V_signal_max * (effective ** self.n_signal) / (
            self.Kd_signal ** self.n_signal + effective ** self.n_signal)
        auto = self.V_auto * (KLF4 ** self.m_klf4) / (
            self.Kd_klf4_auto ** self.m_klf4 + KLF4 ** self.m_klf4)

        dKLF4 = self.V_basal_klf4 + drive + auto - self.k_deg_klf4 * KLF4
        dKRT13 = self.k_krt13 * KLF4 - self.k_deg_krt13 * KRT13
        return [dKLF4, dKRT13]


# ============================================================
# 模型 B: 扩展模型 — 脱敏 + Ca²⁺ + HIF 正反馈
# 状态变量: [KLF4, KRT13, Desens, Ca, HIF]
# ============================================================
class ExtendedModel:
    """
    新增机制:
      Desens: TRPM5 脱敏状态 (0=活性, 1=脱敏)
      Ca:     细胞内钙浓度 (μM), 由 V_m → VGCC 间接调控
      HIF:    缺氧诱导因子, 形成 KLF4↑→化生→缺氧→HIF→KLF4↑ 恶性循环
    """

    def __init__(self):
        # KLF4 核心
        self.k_deg_klf4 = P["k_deg_klf4"]
        self.V_basal_klf4 = P["V_basal_klf4"]
        self.V_auto = P["V_auto"]
        self.m_klf4 = P["m_klf4"]
        self.Kd_klf4_auto = P["Kd_klf4_auto"]

        # TRPM5 信号
        self.V_signal_max = P["V_signal_max"]
        self.n_signal = P["n_signal"]
        self.Kd_signal = P["Kd_signal"]

        # 脱敏模块
        self.k_desens = P["k_desens"] * 60.0   # min⁻¹ → h⁻¹
        self.k_recovery = P["k_recovery"] * 60.0

        # Ca²⁺ 模块
        self.k_VGCC = P["k_VGCC"]
        self.Ca_basal = P["Ca_basal"]
        self.k_Ca_clear = P["k_Ca_clear"]

        # HIF 正反馈
        self.k_hif_klf4 = P["k_hif_klf4"]
        self.Kd_hif = P["Kd_hif"]
        self.k_hif_prod = P["k_hif_prod"]
        self.k_deg_hif = P["k_deg_hif"]

        # KRT13
        self.k_krt13 = P["k_krt13"]
        self.k_deg_krt13 = P["k_deg_krt13"]

    def __call__(self, t, y, TRPM5_input, stress):
        KLF4, KRT13, Desens, Ca, HIF = y

        # --- 1. TRPM5 有效活性 (受脱敏调制) ---
        TRPM5_eff = TRPM5_input * (1.0 - Desens)

        # --- 2. TRPM5 → Ca²⁺ 信号 (Ca 直接受 TRPM5 驱动的 Na⁺ 去极化影响) ---
        # 简化: Na⁺ 去极化程度 ∝ TRPM5_eff, 通过 VGCC 引起 Ca²⁺ 内流
        Ca_influx = self.k_VGCC * TRPM5_eff * stress
        dCa = Ca_influx - self.k_Ca_clear * (Ca - self.Ca_basal)

        # --- 3. TRPM5 脱敏动力学 ---
        # 高 Ca²⁺ 促进脱敏, 低 Ca²⁺ 时恢复
        Ca_norm = Ca / (Ca + 10.0)  # 归一化 Ca 对脱敏的驱动
        dDesens = self.k_desens * TRPM5_eff * Ca_norm - self.k_recovery * Desens

        # --- 4. KLF4 动力学 (多重输入) ---
        # 4a. TRPM5-Ca²⁺ 信号驱动
        signal_drive = self.V_signal_max * (Ca ** self.n_signal) / (
            self.Kd_signal ** self.n_signal + Ca ** self.n_signal)

        # 4b. HIF 正反馈驱动
        hif_drive = self.k_hif_klf4 * (HIF ** 2) / (
            self.Kd_hif ** 2 + HIF ** 2)

        # 4c. KLF4 自维持
        auto = self.V_auto * (KLF4 ** self.m_klf4) / (
            self.Kd_klf4_auto ** self.m_klf4 + KLF4 ** self.m_klf4)

        dKLF4 = (self.V_basal_klf4 + signal_drive + hif_drive + auto
                 - self.k_deg_klf4 * KLF4)

        # --- 5. HIF 动力学 (KLF4→化生→缺氧→HIF) ---
        # 缺氧程度与 KLF4 水平成正比 (更多化生 → 更缺氧)
        hypoxia_signal = KLF4 / (KLF4 + 5.0)
        dHIF = (self.k_hif_prod * hypoxia_signal * stress
                - self.k_deg_hif * HIF)

        # --- 6. KRT13 下游 ---
        dKRT13 = self.k_krt13 * KLF4 - self.k_deg_krt13 * KRT13

        return [dKLF4, dKRT13, dDesens, dCa, dHIF]


# ============================================================
# 仿真运行
# ============================================================
def run_simulation(model, t_span=(0, 80), n_points=1000,
                   TRPM5_levels=None, stress=1.2, y0=None):
    """运行多组 TRPM5 干预水平的仿真, 返回结果字典。"""
    if TRPM5_levels is None:
        TRPM5_levels = [1.0, 0.6, 0.1]
    t_eval = np.linspace(*t_span, n_points)

    results = {}
    for level in TRPM5_levels:
        sol = solve_ivp(model, t_span, y0, t_eval=t_eval,
                        args=(level, stress), method='LSODA')
        results[level] = sol
    return results, t_eval


# ============================================================
# 参数敏感性分析
# ============================================================
def sensitivity_scan(param_name, param_range, model_class=BaseModel,
                     TRPM5_activity=1.0, stress=1.2, y0=None):
    """扫描单个参数对 KLF4 稳态值的影响。"""
    if y0 is None:
        y0 = [5.0, 20.0]
    values = np.linspace(*param_range, 50)
    steady_states = []

    for val in values:
        m = model_class()
        setattr(m, param_name, val)
        sol = solve_ivp(m, (0, 100), y0, t_eval=[100],
                        args=(TRPM5_activity, stress))
        steady_states.append(sol.y[0, -1])

    return values, np.array(steady_states)


def sensitivity_heatmap(param_x, param_y, range_x, range_y,
                        model_class=BaseModel, TRPM5_activity=1.0):
    """二维参数扫描, 返回稳态 KLF4 矩阵用于热力图。"""
    vals_x = np.linspace(*range_x, 30)
    vals_y = np.linspace(*range_y, 30)
    matrix = np.zeros((len(vals_y), len(vals_x)))

    for i, py in enumerate(vals_y):
        for j, px in enumerate(vals_x):
            m = model_class()
            setattr(m, param_x, px)
            setattr(m, param_y, py)
            sol = solve_ivp(m, (0, 80), [5.0, 20.0], t_eval=[80],
                            args=(TRPM5_activity, 1.2))
            matrix[i, j] = sol.y[0, -1]

    return vals_x, vals_y, matrix


# ============================================================
# 主程序
# ============================================================
if __name__ == "__main__":
    stress = 1.2
    t_span = (0, 80)

    # --- 基础模型仿真 ---
    base = BaseModel()
    y0_base = [5.0, 20.0]
    results_base, t_eval = run_simulation(
        base, t_span=t_span, y0=y0_base, stress=stress,
        TRPM5_levels=[1.0, 0.6, 0.1])

    # --- 扩展模型仿真 ---
    ext = ExtendedModel()
    y0_ext = [5.0, 20.0, 0.1, 0.5, 0.2]  # KLF4, KRT13, Desens, Ca, HIF
    results_ext, _ = run_simulation(
        ext, t_span=t_span, y0=y0_ext, stress=stress,
        TRPM5_levels=[1.0, 0.6, 0.1])

    # --- 敏感性分析 ---
    sens_params = [
        ("V_auto", PARAM_RANGES["V_auto"], "KLF4 自维持强度"),
        ("m_klf4", PARAM_RANGES["m_klf4"], "协同系数 m"),
        ("k_deg_klf4", PARAM_RANGES["k_deg_klf4"], "KLF4 降解速率"),
        ("V_signal_max", PARAM_RANGES["V_signal_max"], "信号驱动强度"),
    ]

    # ============================================================
    # 可视化
    # ============================================================
    fig = plt.figure(figsize=(18, 14))

    # --- Row 1: 基础模型 KLF4 + KRT13 ---
    ax1 = fig.add_subplot(3, 4, 1)
    colors = {1.0: '#d62728', 0.6: '#ff7f0e', 0.1: '#2ca02c'}
    labels = {1.0: 'TRPM5 100% (无干预)', 0.6: 'TRPM5 60% (弱干预)',
              0.1: 'TRPM5 10% (强阻断)'}
    for level, sol in results_base.items():
        ax1.plot(sol.t, sol.y[0], color=colors[level], label=labels[level])
    ax1.axhspan(0, 2.0, color='green', alpha=0.08)
    ax1.set_title('Base: KLF4 Dynamics')
    ax1.set_ylabel('KLF4')
    ax1.legend(fontsize=7)
    ax1.grid(alpha=0.3)

    ax2 = fig.add_subplot(3, 4, 2)
    for level, sol in results_base.items():
        ax2.plot(sol.t, sol.y[1], color=colors[level])
    ax2.set_title('Base: KRT13 Trajectory')
    ax2.set_ylabel('KRT13')
    ax2.set_xlabel('Time (h)')
    ax2.grid(alpha=0.3)

    # --- Row 1 continued: 扩展模型 KLF4 + KRT13 ---
    ax3 = fig.add_subplot(3, 4, 3)
    for level, sol in results_ext.items():
        ax3.plot(sol.t, sol.y[0], color=colors[level])
    ax3.axhspan(0, 2.0, color='green', alpha=0.08)
    ax3.set_title('Extended: KLF4 Dynamics')
    ax3.set_xlabel('Time (h)')
    ax3.grid(alpha=0.3)

    ax4 = fig.add_subplot(3, 4, 4)
    for level, sol in results_ext.items():
        ax4.plot(sol.t, sol.y[1], color=colors[level])
    ax4.set_title('Extended: KRT13 Trajectory')
    ax4.set_xlabel('Time (h)')
    ax4.grid(alpha=0.3)

    # --- Row 2: 扩展模型新增状态变量 ---
    ax5 = fig.add_subplot(3, 4, 5)
    for level, sol in results_ext.items():
        ax5.plot(sol.t, sol.y[2], color=colors[level])
    ax5.set_title('TRPM5 Desensitization')
    ax5.set_ylabel('Desens (0=active)')
    ax5.grid(alpha=0.3)

    ax6 = fig.add_subplot(3, 4, 6)
    for level, sol in results_ext.items():
        ax6.plot(sol.t, sol.y[3], color=colors[level])
    ax6.set_title('[Ca²⁺]ᵢ Dynamics')
    ax6.set_ylabel('Ca²⁺ (μM)')
    ax6.grid(alpha=0.3)

    ax7 = fig.add_subplot(3, 4, 7)
    for level, sol in results_ext.items():
        ax7.plot(sol.t, sol.y[4], color=colors[level])
    ax7.set_title('HIF Positive Feedback')
    ax7.set_ylabel('HIF')
    ax7.set_xlabel('Time (h)')
    ax7.grid(alpha=0.3)

    # 脱敏 vs KLF4 相图
    ax8 = fig.add_subplot(3, 4, 8)
    sol_ref = results_ext[1.0]
    ax8.plot(sol_ref.y[2], sol_ref.y[0], color='#9467bd')
    ax8.set_title('Desens–KLF4 Phase Portrait')
    ax8.set_xlabel('Desensitization')
    ax8.set_ylabel('KLF4')
    ax8.grid(alpha=0.3)

    # --- Row 3: 参数敏感性分析 (4个关键参数) ---
    for idx, (pname, prange, plabel) in enumerate(sens_params):
        ax = fig.add_subplot(3, 4, 9 + idx)
        vals, steady = sensitivity_scan(pname, prange)
        ax.plot(vals, steady, color='#1f77b4', linewidth=2)
        ax.axhline(y=2.0, color='green', linestyle=':', alpha=0.5,
                   label='健康阈值')
        ax.set_title(f'Sensitivity: {plabel}')
        ax.set_xlabel(pname)
        ax.set_ylabel('KLF4 steady state')
        ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.show()

    # --- 打印关键结果 ---
    print("=" * 60)
    print("TRPM5-KLF4 模型 v2.0 — 分析报告")
    print("=" * 60)

    # 基础模型稳态
    print("\n[基础模型] KLF4 稳态值:")
    for level, sol in results_base.items():
        print(f"  TRPM5={level:.0%}: KLF4={sol.y[0,-1]:.2f}, KRT13={sol.y[1,-1]:.2f}")

    # 扩展模型稳态
    print("\n[扩展模型] 稳态值:")
    for level, sol in results_ext.items():
        print(f"  TRPM5={level:.0%}: KLF4={sol.y[0,-1]:.2f}, KRT13={sol.y[1,-1]:.2f}, "
              f"Desens={sol.y[2,-1]:.2f}, Ca={sol.y[3,-1]:.2f}, HIF={sol.y[4,-1]:.2f}")

    # 逆转时间估算
    print("\n[治疗意义] KLF4 回落至 <2.0 所需时间:")
    for level, sol in results_ext.items():
        if level < 0.5:
            mask = sol.y[0] < 2.0
            if mask.any():
                t_reversal = sol.t[mask][0]
                print(f"  TRPM5={level:.0%}: ~{t_reversal:.1f} h")
    print("\n模型采用参数来源于报告第二节真实实验数据。")
