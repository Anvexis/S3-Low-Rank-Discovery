# plots/visualize_ranks.py
"""
Generates publication-ready plots from experiment CSVs.
Requires: pandas, matplotlib
Run from repo root: python3 plots/visualize_ranks.py
"""
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import glob

plt.rcParams.update({
    "font.size": 11,
    "font.family": "serif",
    "axes.grid": True,
    "grid.alpha": 0.3,
    "savefig.bbox": "tight",
    "savefig.dpi": 300
})

RESULTS_DIR = Path("results")
PLOTS_DIR = Path("plots")
PLOTS_DIR.mkdir(exist_ok=True)

def load_csv_safe(pattern):
    files = list(RESULTS_DIR.glob(pattern))
    if not files: return pd.DataFrame()
    return pd.concat([pd.read_csv(f) for f in files], ignore_index=True)

def plot_rank_stability():
    df = load_csv_safe("verified_streaming_rank.csv")
    if df.empty: print("⚠️ No streaming rank data found."); return
    
    plt.figure(figsize=(6, 4))
    plt.plot(df["n"], df["rank"], "o-", color="#2563eb", linewidth=2, markersize=6)
    plt.axhline(y=3, color="#dc2626", linestyle="--", label="Theoretical stable rank = 3")
    plt.xlabel(r"Factor Base Size $|\mathcal{B}|$")
    plt.ylabel("Matrix Rank")
    plt.title(r"Rank Stability of $S_3(x_i, x_j, x_R)$ over $\mathbb{F}_p$")
    plt.legend()
    plt.savefig(PLOTS_DIR / "rank_stability.pdf")
    plt.close()
    print("✅ Saved: plots/rank_stability.pdf")

def plot_scaling():
    df = load_csv_safe("verified_streaming_rank.csv")
    if df.empty: print("⚠️ No scaling data found."); return
    
    fig, ax1 = plt.subplots(figsize=(6, 4))
    color_t = "#059669"
    ax1.set_xlabel(r"Factor Base Size $|\mathcal{B}|$")
    ax1.set_ylabel("Time (s)", color=color_t)
    ax1.plot(df["n"], df["time_s"], "s-", color=color_t, linewidth=2, markersize=6)
    ax1.tick_params(axis='y', labelcolor=color_t)
    
    ax2 = ax1.twinx()
    color_m = "#7c3aed"
    ax2.set_ylabel("Δ Memory (MB)", color=color_m)
    ax2.plot(df["n"], df["mem_delta_MB"], "D--", color=color_m, linewidth=2, markersize=6)
    ax2.tick_params(axis='y', labelcolor=color_m)
    
    plt.title("Streaming Solver Scaling (p ≈ 2⁶⁴)")
    plt.savefig(PLOTS_DIR / "scaling_performance.pdf")
    plt.close()
    print("✅ Saved: plots/scaling_performance.pdf")

def plot_poly_ranks():
    s5 = load_csv_safe("S5_rank_test.csv")
    s6 = load_csv_safe("S6_rank_test.csv")
    if s5.empty and s6.empty: print("⚠️ No S5/S6 data found."); return
    
    plt.figure(figsize=(6, 4))
    if not s5.empty:
        s5_data = s5[s5["poly"] == "S5"]
        plt.plot(s5_data["n"], s5_data["rank"], "o-", label="$S_5$ (fixed 3 targets)", color="#ea580c")
    if not s6.empty:
        s6_data = s6[s6["poly"] == "S6"]
        plt.plot(s6_data["n"], s6_data["rank"], "s-", label="$S_6$ (fixed 4 targets)", color="#9333ea")
        
    plt.axhline(y=3, color="#2563eb", linestyle=":", label="$S_3$ baseline = 3")
    plt.xlabel(r"Factor Base Size $|\mathcal{B}|$")
    plt.ylabel("Observed Rank")
    plt.title("Higher Summation Polynomials: Rank Growth")
    plt.legend()
    plt.savefig(PLOTS_DIR / "poly_rank_comparison.pdf")
    plt.close()
    print("✅ Saved: plots/poly_rank_comparison.pdf")

def main():
    print("📊 Generating publication plots...")
    plot_rank_stability()
    plot_scaling()
    plot_poly_ranks()
    print("\n🎨 All plots saved to plots/")

if __name__ == "__main__":
    main()