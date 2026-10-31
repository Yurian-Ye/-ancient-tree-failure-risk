# -*- coding: utf-8 -*-
"""
Wind load vs root anchorage: safety factor of ancient trees
===========================================================
A minimal, reproducible model used to check whether the "structure -> safety
factor -> risk class" chain in our proposal actually runs, and to decide which
field parameters should be measured first.

Model
-----
    q        = 0.5 * rho * v^2                 wind pressure          [Pa]
    A        = 0.65 * (0.75 H) * W             effective crown area   [m2]
    M_wind   = Cd * q * A * (0.6 H) * (1+lean) overturning moment     [kN m]
    M_root   = 0.3 * su * D * h_eff^2          root-plate resistance   [kN m]
               h_eff = max(h - d, 0)
    SF       = M_root / M_wind

All parameters are order-of-magnitude values taken from public design codes and
the literature. **None of them is calibrated against field measurements.**
The model is meant for ranking and sensitivity analysis, not for judging the
safety of an individual tree.

Outputs (written to ./figures)
------------------------------
    fig1_sf_cdf.png            cumulative distribution of SF (Monte Carlo)
    fig2_scour_vs_sf.png       SF versus scour depth around the root plate
    fig3_sobol.png             Sobol first-order / total-order indices
    metrics.json               every number quoted in the README

Run
---
    pip install -r requirements.txt
    python src/mech_model.py
"""
import json
import os

import numpy as np
from scipy.stats import gumbel_r

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "figures")
os.makedirs(OUT, exist_ok=True)

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

# ----------------------------------------------------------------- reference tree
H_BASE = 22.0      # tree height [m]
W_BASE = 17.0      # crown width [m]
LEAN_BASE = 0.15   # crown eccentricity / crown width

RHO = 1.225        # air density [kg/m3]
W0 = 0.70          # 50-year basic wind pressure, Fuzhou [kN/m2] (GB 50009 order)
V50 = (2 * W0 * 1000 / RHO) ** 0.5   # 10-min mean wind speed, ~33.8 m/s
CV = 0.15                            # coefficient of variation of annual maxima
SIG = CV * V50
BETA = SIG * (6 ** 0.5) / np.pi      # Gumbel scale
MU = V50 - BETA * (-np.log(-np.log(1 - 1 / 50)))   # Gumbel location

PARAMS = {
    "v":  dict(lo=12.0, hi=60.0, label="极值风速"),
    "Cd": dict(lo=0.15, hi=0.30, label="冠层阻力系数"),
    "D":  dict(lo=8.0,  hi=16.0, label="根盘直径"),
    "h":  dict(lo=1.5,  hi=3.5,  label="根盘埋深"),
    "su": dict(lo=25.0, hi=80.0, label="土体抗剪强度"),
    "d":  dict(lo=0.0,  hi=3.0,  label="树基冲刷深度"),
}
NAMES = list(PARAMS)


def sample(name, n, rng):
    """Draw n samples of one parameter from its marginal distribution."""
    p = PARAMS[name]
    if name == "v":
        out = np.empty(n)
        m = 0
        while m < n:
            cand = gumbel_r.rvs(loc=MU, scale=BETA, size=n - m, random_state=rng)
            cand = cand[(cand >= p["lo"]) & (cand <= p["hi"])]
            out[m:m + len(cand)] = cand
            m += len(cand)
        return out
    if name == "d":
        # Most trees show no scour; a minority on river banks and slope toes do.
        u = rng.random(n)
        shallow = rng.uniform(0.0, 0.5, n)
        deep = rng.uniform(0.5, 3.0, n)
        return np.where(u < 0.8, shallow, deep)
    return rng.uniform(p["lo"], p["hi"], n)


def model(v, Cd, D, h, su, d, H=H_BASE, W=W_BASE, lean=LEAN_BASE):
    """Safety factor = root anchorage moment / wind overturning moment."""
    q = 0.5 * RHO * v ** 2 / 1000.0          # [kN/m2]
    A = 0.65 * 0.75 * H * W                   # [m2]
    M_wind = Cd * q * A * (0.6 * H) * (1 + lean)
    h_eff = np.clip(h - d, 0.0, None)
    M_root = 0.3 * su * D * h_eff ** 2
    return M_root / np.maximum(M_wind, 1e-9)


def monte_carlo(n=20000, seed=42):
    rng = np.random.default_rng(seed)
    cols = [sample(k, n, rng) for k in NAMES]
    return model(*cols), dict(zip(NAMES, cols))


def sobol(n=2048, seed=7):
    """Saltelli (2010) estimators; inputs follow their own marginals."""
    rng = np.random.default_rng(seed)
    d = len(NAMES)
    A = np.column_stack([sample(k, n, rng) for k in NAMES])
    B = np.column_stack([sample(k, n, rng) for k in NAMES])
    fA, fB = model(*A.T), model(*B.T)
    var = np.var(np.concatenate([fA, fB]))
    s1, st = [], []
    for i in range(d):
        AB = A.copy()
        AB[:, i] = B[:, i]
        fAB = model(*AB.T)
        s1.append(np.mean(fB * (fAB - fA)) / var)
        st.append(0.5 * np.mean((fA - fAB) ** 2) / var)
    return np.array(s1), np.array(st), n * (d + 2)


def main():
    sf, draws = monte_carlo()
    s1, st, n_runs = sobol()

    med = {k: float(np.median(draws[k])) for k in NAMES}
    sf_bench = float(model(**med))
    sf_design = float(model(v=V50, Cd=0.225, D=12.0, h=2.5, su=52.5, d=0.0))

    lo, hi = 0.0, 3.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if model(**{**med, "d": mid}) > 1.0:
            lo = mid
        else:
            hi = mid
    d_crit = 0.5 * (lo + hi)

    share_high = float(np.mean(sf < 1.0))
    share_mid = float(np.mean((sf >= 1.0) & (sf < 1.5)))

    res = dict(
        v50=round(V50, 2), n_mc=int(sf.size), n_sobol_runs=int(n_runs),
        sf_median_bench=round(sf_bench, 3), sf_design_wind=round(sf_design, 3),
        sf_p05=round(float(np.percentile(sf, 5)), 3),
        sf_p50=round(float(np.percentile(sf, 50)), 3),
        sf_p95=round(float(np.percentile(sf, 95)), 3),
        share_high_risk=round(share_high, 4),
        share_mid_risk=round(share_mid, 4),
        d_crit=round(d_crit, 3),
        sobol={k: dict(S1=round(float(a), 3), ST=round(float(b), 3))
               for k, a, b in zip(NAMES, s1, st)},
        base=dict(H=H_BASE, W=W_BASE, lean=LEAN_BASE, W0=W0),
    )
    with open(os.path.join(OUT, "metrics.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)

    # ---- fig 1: cumulative distribution
    nz = sf[sf > 0]
    share0 = 1 - nz.size / sf.size
    xs = np.sort(nz)
    ys = np.arange(1, xs.size + 1) / xs.size

    fig, ax = plt.subplots(figsize=(7.2, 4.2), dpi=200)
    ax.semilogx(xs, ys, lw=2.0, color="#2f4858")
    ax.axvline(1.0, color="#b23a2e", lw=1.5, ls="--")
    ax.axvline(1.5, color="#c8862a", lw=1.5, ls="--")
    ax.set_xlim(0.02, 30)
    ax.set_ylim(0, 1)
    ax.set_xlabel("安全系数 SF（对数刻度，根系锚固抗力 / 风致倾覆弯矩）")
    ax.set_ylabel("累积占比")
    ax.set_title(f"蒙特卡洛传播：{sf.size:,} 组工况的安全系数累积分布", fontsize=11, pad=10)
    ax.grid(alpha=0.25, lw=0.6)
    ax.annotate(f"高风险 {share_high:.1%}", xy=(1.0, share_high), xytext=(0.055, 0.70),
                fontsize=9, color="#b23a2e",
                arrowprops=dict(arrowstyle="->", color="#b23a2e", lw=1.0))
    ax.annotate(f"中风险及以下 {(share_high + share_mid):.1%}",
                xy=(1.5, share_high + share_mid), xytext=(3.1, 0.60),
                fontsize=9, color="#c8862a",
                arrowprops=dict(arrowstyle="->", color="#c8862a", lw=1.0))
    ax.text(0.035, 0.94, f"另有 {share0:.1%} 工况根盘完全悬空\n（冲刷深度超过埋深，SF 记为 0）",
            transform=ax.transAxes, ha="left", va="top", fontsize=8.2, color="#7a6c53")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig1_sf_cdf.png"))
    plt.close(fig)

    # ---- fig 2: scour depth
    ds = np.linspace(0, 3, 121)
    fig, ax = plt.subplots(figsize=(7.2, 4.0), dpi=200)
    for lean, ls in [(0.0, ":"), (LEAN_BASE, "-"), (0.35, "--")]:
        ys2 = [model(**{**med, "d": x}, lean=lean) for x in ds]
        ax.plot(ds, ys2, lw=1.9, ls=ls, color="#2f4858", label=f"偏冠率 {lean:.0%}")
    ax.axhline(1.5, color="#c8862a", lw=1.4, ls="--")
    ax.axhline(1.0, color="#b23a2e", lw=1.4, ls="--")
    ax.axvline(d_crit, color="#7a6c53", lw=1.2, ls="-.")
    ax.text(d_crit + 0.06, 5.5, f"临界冲刷深度\n≈ {d_crit:.2f} m",
            color="#7a6c53", fontsize=9, va="top")
    ax.set_xlabel("树基冲刷深度 d（m）")
    ax.set_ylabel("安全系数 SF")
    ax.set_ylim(0, 7.0)
    ax.set_title("其余参数取中位时，安全系数随树基冲刷深度的变化", fontsize=11, pad=10)
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig2_scour_vs_sf.png"))
    plt.close(fig)

    # ---- fig 3: Sobol
    order = np.argsort(-st)
    labels = [PARAMS[NAMES[i]]["label"] for i in order]
    fig, ax = plt.subplots(figsize=(7.2, 3.6), dpi=200)
    y = np.arange(len(order))
    ax.barh(y - 0.19, s1[order], height=0.36, color="#7ba189", label="一阶 S1")
    ax.barh(y + 0.19, st[order], height=0.36, color="#2f4858", label="总效应 ST")
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("Sobol 敏感性指数")
    ax.set_title(f"Sobol 全局敏感性（Saltelli 估计量，{n_runs:,} 次模型调用）", fontsize=11, pad=10)
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    for i, v in enumerate(st[order]):
        ax.text(v + 0.012, i + 0.19, f"{v:.2f}", va="center", fontsize=8.5, color="#2f4858")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig3_sobol.png"))
    plt.close(fig)

    print(json.dumps({k: v for k, v in res.items() if k != "base"},
                     ensure_ascii=False, indent=1))
    print("figures ->", OUT)


if __name__ == "__main__":
    main()
