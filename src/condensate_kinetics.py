#!/usr/bin/env python3
"""Reaction kinetics inside biomolecular condensates.

Models the concentration enhancement effect — the key argument for
condensates as protocellular compartments in origin-of-life scenarios.

Inside a condensate, molecules partition preferentially:
    c_in = K · c_out     (K = partition coefficient, typically 10–1000)

For a bimolecular reaction A + B → C with rate constant k:
    rate_dilute    = k · [A]_out · [B]_out
    rate_condensate = k · K_A · [A]_out · K_B · [B]_out = K_A · K_B · rate_dilute

The reaction rate is enhanced by K² (if K_A ≈ K_B ≈ K), which can be
10⁴–10⁶ fold. This is how condensates could have kickstarted prebiotic
chemistry in an otherwise impossibly dilute primordial ocean.

We also model Michaelis-Menten ribozyme kinetics inside coacervate
droplets, showing conditions where replication exceeds degradation
only within the condensate.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

FIGURES = Path(__file__).resolve().parent.parent / "figures"
FIGURES.mkdir(exist_ok=True)


def partition_enhancement(K_A, K_B=None):
    """Rate enhancement factor for A + B → C inside a condensate."""
    if K_B is None:
        K_B = K_A
    return K_A * K_B


def michaelis_menten(S, Vmax, Km):
    """Michaelis-Menten reaction rate."""
    return Vmax * S / (Km + S)


def ribozyme_replication_rate(S, K, Vmax=1.0, Km=10.0):
    """Ribozyme-catalysed replication rate inside a condensate.

    The substrate concentration inside the condensate is K·S,
    and the enzyme concentration is also enhanced by K, giving
    effective Vmax_eff = K · Vmax.
    """
    S_eff = K * S
    Vmax_eff = K * Vmax
    return michaelis_menten(S_eff, Vmax_eff, Km)


def degradation_rate(S, K, k_deg=0.1):
    """First-order degradation rate inside condensate."""
    return k_deg * K * S


def plot_rate_enhancement():
    """Show how reaction rate scales with partition coefficient."""
    K = np.logspace(0, 3, 200)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

    # Bimolecular enhancement
    ax1.loglog(K, K**2, 'b-', lw=2.5, label='Bimolecular (A + B → C)')
    ax1.loglog(K, K, 'g--', lw=2, label='Unimolecular (A → B)')
    ax1.loglog(K, K**3, 'r:', lw=2, label='Trimolecular (A + B + C → D)')

    ax1.axhspan(1e4, 1e6, alpha=0.1, color='orange',
                label='Biologically relevant range')
    ax1.axvline(100, color='grey', ls=':', lw=1, alpha=0.5)
    ax1.annotate('Typical condensate\nK ≈ 100', xy=(100, 1),
                 fontsize=9, ha='center', va='bottom',
                 xytext=(100, 3), arrowprops=dict(arrowstyle='->', color='grey'))

    ax1.set_xlabel("Partition coefficient K", fontsize=13)
    ax1.set_ylabel("Rate enhancement factor", fontsize=13)
    ax1.set_title("Reaction rate enhancement in condensates", fontsize=14)
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3, which='both')
    ax1.set_xlim(1, 1000)
    ax1.set_ylim(1, 1e9)

    # Effective concentration
    c_out = np.logspace(-9, -3, 200)  # molar
    for Ki in [1, 10, 100, 1000]:
        ax2.loglog(c_out * 1e6, Ki * c_out * 1e6, lw=2, label=f'K = {Ki}')

    ax2.axhline(1e3, color='grey', ls=':', lw=1)
    ax2.text(1e-1, 1.5e3, 'mM range (enzyme optimal)', fontsize=9, color='grey')
    ax2.axhline(1, color='grey', ls=':', lw=1)
    ax2.text(1e-1, 1.5, 'μM range (typical Km)', fontsize=9, color='grey')

    ax2.set_xlabel("Bulk concentration (μM)", fontsize=13)
    ax2.set_ylabel("Effective concentration inside condensate (μM)", fontsize=13)
    ax2.set_title("Concentration enhancement", fontsize=14)
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3, which='both')

    fig.tight_layout()
    fig.savefig(FIGURES / "rate_enhancement.png", dpi=200)
    plt.close(fig)
    print(f"  saved {FIGURES / 'rate_enhancement.png'}")


def plot_ribozyme_kinetics():
    """Show Michaelis-Menten kinetics inside vs outside condensate."""
    S = np.linspace(0, 50, 300)  # substrate concentration (μM)
    Vmax = 1.0  # μM/min
    Km = 10.0   # μM

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

    # Panel 1: reaction rates at different K
    for K in [1, 5, 20, 100]:
        v = ribozyme_replication_rate(S, K, Vmax, Km)
        ax1.plot(S, v, lw=2, label=f'K = {K}')

    ax1.set_xlabel("[Substrate]_out  (μM)", fontsize=13)
    ax1.set_ylabel("Replication rate  (μM/min)", fontsize=13)
    ax1.set_title("Ribozyme replication rate inside condensate", fontsize=14)
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)

    # Panel 2: replication vs degradation — the origin-of-life window
    K_values = np.linspace(1, 200, 300)
    S_fixed = 1.0  # 1 μM bulk substrate
    k_deg = 0.05

    v_rep = np.array([ribozyme_replication_rate(S_fixed, K, Vmax, Km) for K in K_values])
    v_deg = np.array([degradation_rate(S_fixed, K, k_deg) for K in K_values])

    ax2.plot(K_values, v_rep, 'b-', lw=2.5, label='Replication rate')
    ax2.plot(K_values, v_deg, 'r--', lw=2.5, label='Degradation rate')
    ax2.fill_between(K_values, v_rep, v_deg,
                     where=v_rep > v_deg, alpha=0.15, color='green',
                     label='Net growth window')
    ax2.fill_between(K_values, v_rep, v_deg,
                     where=v_rep <= v_deg, alpha=0.1, color='red',
                     label='Net decay')

    # Find crossover
    crossover_idx = np.where(np.diff(np.sign(v_rep - v_deg)))[0]
    if len(crossover_idx) > 0:
        K_cross = K_values[crossover_idx[0]]
        ax2.axvline(K_cross, color='green', ls=':', lw=1.5)
        ax2.annotate(f'K* = {K_cross:.0f}\n(critical partition\n coefficient)',
                     xy=(K_cross, v_rep[crossover_idx[0]]),
                     xytext=(K_cross + 30, v_rep[crossover_idx[0]] * 1.5),
                     fontsize=10, arrowprops=dict(arrowstyle='->', color='green'))

    ax2.set_xlabel("Partition coefficient K", fontsize=13)
    ax2.set_ylabel("Rate  (μM/min)", fontsize=13)
    ax2.set_title("Replication vs degradation: origin-of-life window", fontsize=14)
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(FIGURES / "ribozyme_kinetics.png", dpi=200)
    plt.close(fig)
    print(f"  saved {FIGURES / 'ribozyme_kinetics.png'}")


def plot_dilute_vs_condensate():
    """Side-by-side comparison of reaction landscape."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    c_bulk = np.logspace(-3, 1, 200)  # μM
    K = 100
    k2 = 0.01  # bimolecular rate constant

    # Panel 1: dilute solution
    rate_dilute = k2 * c_bulk**2
    axes[0].loglog(c_bulk, rate_dilute, 'b-', lw=2.5)
    axes[0].set_title("Dilute solution", fontsize=14, fontweight='bold')
    axes[0].set_xlabel("[Reactant] (μM)", fontsize=12)
    axes[0].set_ylabel("Reaction rate (μM/s)", fontsize=12)
    axes[0].axhline(1e-2, color='orange', ls='--', lw=1.5, label='Biological threshold')
    axes[0].axvspan(0.001, 0.1, alpha=0.1, color='blue', label='Prebiotic ocean\n([X] ~ 1–100 nM)')
    axes[0].legend(fontsize=9)
    axes[0].grid(True, alpha=0.3, which='both')
    axes[0].set_ylim(1e-10, 1e2)

    # Panel 2: inside condensate
    rate_cond = k2 * (K * c_bulk)**2
    axes[1].loglog(c_bulk, rate_cond, 'r-', lw=2.5)
    axes[1].set_title(f"Inside condensate (K = {K})", fontsize=14, fontweight='bold')
    axes[1].set_xlabel("[Reactant]_bulk (μM)", fontsize=12)
    axes[1].set_ylabel("Reaction rate (μM/s)", fontsize=12)
    axes[1].axhline(1e-2, color='orange', ls='--', lw=1.5, label='Biological threshold')
    axes[1].axvspan(0.001, 0.1, alpha=0.1, color='blue', label='Prebiotic ocean')
    axes[1].legend(fontsize=9)
    axes[1].grid(True, alpha=0.3, which='both')
    axes[1].set_ylim(1e-10, 1e2)

    # Panel 3: enhancement factor
    enhancement = rate_cond / rate_dilute
    axes[2].loglog(c_bulk, enhancement, 'purple', lw=2.5)
    axes[2].set_ylim(K**2 * 0.5, K**2 * 2)
    axes[2].axhline(K**2, color='grey', ls=':', lw=1.5)
    axes[2].text(0.01, K**2 * 1.15, f'K² = {K**2:.0e}', fontsize=11, color='grey')
    axes[2].set_title("Enhancement factor", fontsize=14, fontweight='bold')
    axes[2].set_xlabel("[Reactant]_bulk (μM)", fontsize=12)
    axes[2].set_ylabel("Rate_condensate / Rate_dilute", fontsize=12)
    axes[2].grid(True, alpha=0.3, which='both')

    fig.suptitle("Condensate-assisted chemistry: the concentration advantage",
                 fontsize=15, y=1.02)
    fig.tight_layout()
    fig.savefig(FIGURES / "dilute_vs_condensate.png", dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"  saved {FIGURES / 'dilute_vs_condensate.png'}")


def plot_prebiotic_phase_space():
    """Phase diagram: which combinations of K and bulk concentration
    allow net replication."""
    K_range = np.logspace(0, 3, 200)
    c_range = np.logspace(-3, 1, 200)
    KK, CC = np.meshgrid(K_range, c_range)

    Vmax = 1.0
    Km = 10.0
    k_deg = 0.05

    v_rep = KK * Vmax * (KK * CC) / (Km + KK * CC)
    v_deg = k_deg * KK * CC
    net = v_rep - v_deg

    fig, ax = plt.subplots(figsize=(8, 6))
    cs = ax.contourf(KK, CC, net, levels=20, cmap='RdBu', extend='both')
    ax.contour(KK, CC, net, levels=[0], colors='black', linewidths=2.5)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel("Partition coefficient K", fontsize=13)
    ax.set_ylabel("Bulk substrate concentration (μM)", fontsize=13)
    ax.set_title("Net replication rate: replication − degradation", fontsize=14)
    cbar = fig.colorbar(cs, label='Net rate (μM/min)')

    ax.text(200, 0.005, 'Net growth\n(life possible)', fontsize=12,
            ha='center', color='navy', fontweight='bold')
    ax.text(3, 5, 'Net decay', fontsize=12,
            ha='center', color='darkred', fontweight='bold')

    fig.tight_layout()
    fig.savefig(FIGURES / "prebiotic_phase_space.png", dpi=200)
    plt.close(fig)
    print(f"  saved {FIGURES / 'prebiotic_phase_space.png'}")


if __name__ == "__main__":
    print("Condensate reaction kinetics")
    print("=" * 40)

    K = 100
    print(f"  Partition coefficient K = {K}")
    print(f"  Bimolecular enhancement: K² = {K**2:.0e}")
    print(f"  At 1 nM bulk: effective concentration = {K * 1e-3:.1f} μM")

    print("\nGenerating figures...")
    plot_rate_enhancement()
    plot_ribozyme_kinetics()
    plot_dilute_vs_condensate()
    plot_prebiotic_phase_space()
    print("Done.")
