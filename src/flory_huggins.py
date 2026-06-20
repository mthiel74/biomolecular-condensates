#!/usr/bin/env python3
"""Flory-Huggins lattice theory of polymer-solvent phase separation.

Implements the free energy functional, spinodal/binodal calculations,
and phase diagram construction for biomolecular condensate formation.

The Flory-Huggins free energy per lattice site:

    f(φ) = φ ln(φ)/N₁ + (1-φ) ln(1-φ)/N₂ + χ φ(1-φ)

where φ is the polymer volume fraction, N₁ and N₂ are the degrees
of polymerisation of each species, and χ is the Flory interaction
parameter (χ ~ 1/T for upper critical solution temperature systems).
"""

import numpy as np
from scipy.optimize import brentq, fsolve
import matplotlib.pyplot as plt
from pathlib import Path

FIGURES = Path(__file__).resolve().parent.parent / "figures"
FIGURES.mkdir(exist_ok=True)


def free_energy_density(phi, chi, N1=100, N2=1):
    """Flory-Huggins free energy density f(φ)/kT per lattice site."""
    phi = np.asarray(phi, dtype=float)
    f = np.full_like(phi, np.inf)
    mask = (phi > 0) & (phi < 1)
    p = phi[mask]
    f[mask] = (p * np.log(p) / N1
               + (1 - p) * np.log(1 - p) / N2
               + chi * p * (1 - p))
    return f


def chemical_potential(phi, chi, N1=100, N2=1):
    """Exchange chemical potential μ = ∂f/∂φ."""
    phi = np.asarray(phi, dtype=float)
    return ((np.log(phi) + 1) / N1
            - (np.log(1 - phi) + 1) / N2
            + chi * (1 - 2 * phi))


def osmotic_pressure(phi, chi, N1=100, N2=1):
    """Osmotic pressure Π = f - φ μ  (used for common tangent)."""
    f = free_energy_density(phi, chi, N1, N2)
    mu = chemical_potential(phi, chi, N1, N2)
    return f - phi * mu


def d2f_dphi2(phi, chi, N1=100, N2=1):
    """Second derivative ∂²f/∂φ² — zero on the spinodal."""
    return 1.0 / (N1 * phi) + 1.0 / (N2 * (1 - phi)) - 2 * chi


def critical_point(N1=100, N2=1):
    """Analytic critical point (φ_c, χ_c)."""
    sqrt_ratio = np.sqrt(N2 / N1)
    phi_c = sqrt_ratio / (1 + sqrt_ratio)
    chi_c = 0.5 * (1.0 / np.sqrt(N1) + 1.0 / np.sqrt(N2))**2
    return phi_c, chi_c


def spinodal_curve(N1=100, N2=1, n_points=200):
    """Spinodal curve: locus of ∂²f/∂φ² = 0."""
    phi_c, chi_c = critical_point(N1, N2)
    eps = 1e-6
    phi_vals = np.linspace(eps, 1 - eps, n_points)
    chi_spinodal = 0.5 * (1.0 / (N1 * phi_vals) + 1.0 / (N2 * (1 - phi_vals)))
    mask = chi_spinodal >= chi_c
    return phi_vals[mask], chi_spinodal[mask]


def binodal_curve(N1=100, N2=1, n_points=80):
    """Binodal curve via common tangent construction.

    At each χ > χ_c, find φ_L and φ_R such that
    μ(φ_L) = μ(φ_R) and Π(φ_L) = Π(φ_R).
    """
    phi_c, chi_c = critical_point(N1, N2)
    chi_max = chi_c * 3.0
    chi_vals = np.linspace(chi_c + 1e-4, chi_max, n_points)

    phi_left = []
    phi_right = []
    chi_good = []

    for chi in chi_vals:
        spinodal_roots = _spinodal_roots(chi, N1, N2)
        if spinodal_roots is None:
            continue
        sp_l, sp_r = spinodal_roots

        def equations(x):
            pL, pR = x
            return [
                chemical_potential(pL, chi, N1, N2) - chemical_potential(pR, chi, N1, N2),
                osmotic_pressure(pL, chi, N1, N2) - osmotic_pressure(pR, chi, N1, N2)
            ]

        try:
            sol = fsolve(equations, [sp_l * 0.5, sp_r + (1 - sp_r) * 0.5],
                         full_output=True)
            x_sol, info, ier, msg = sol
            if ier == 1 and 0 < x_sol[0] < x_sol[1] < 1:
                phi_left.append(x_sol[0])
                phi_right.append(x_sol[1])
                chi_good.append(chi)
        except Exception:
            continue

    return np.array(chi_good), np.array(phi_left), np.array(phi_right)


def _spinodal_roots(chi, N1, N2):
    """Find the two roots of ∂²f/∂φ² = 0 at given χ."""
    phi_c, chi_c = critical_point(N1, N2)
    if chi < chi_c:
        return None
    eps = 1e-8
    try:
        r1 = brentq(lambda p: d2f_dphi2(p, chi, N1, N2), eps, phi_c)
        r2 = brentq(lambda p: d2f_dphi2(p, chi, N1, N2), phi_c, 1 - eps)
        return r1, r2
    except ValueError:
        return None


def plot_free_energy(N1=100, N2=1, chi_values=None):
    """Plot f(φ) for several χ values."""
    if chi_values is None:
        _, chi_c = critical_point(N1, N2)
        chi_values = [0.3 * chi_c, 0.6 * chi_c, chi_c, 1.5 * chi_c, 2.5 * chi_c]

    phi = np.linspace(1e-4, 1 - 1e-4, 500)
    _, chi_c = critical_point(N1, N2)

    fig, ax = plt.subplots(figsize=(8, 5))
    colors = plt.cm.coolwarm(np.linspace(0.1, 0.9, len(chi_values)))

    for chi, c in zip(chi_values, colors):
        f = free_energy_density(phi, chi, N1, N2)
        label = f"χ = {chi:.3f}" + (" (critical)" if abs(chi - chi_c) < 1e-4 else "")
        ax.plot(phi, f, color=c, lw=2, label=label)

    ax.set_xlabel("Volume fraction φ", fontsize=13)
    ax.set_ylabel("Free energy density  f(φ) / kT", fontsize=13)
    ax.set_title(f"Flory–Huggins free energy  (N₁ = {N1}, N₂ = {N2})", fontsize=14)
    ax.legend(fontsize=10, framealpha=0.9)
    ax.set_xlim(0, 1)
    ymin = np.nanmin(free_energy_density(phi, max(chi_values), N1, N2))
    ax.set_ylim(ymin - 0.05, 0.02)
    ax.axhline(0, color='grey', lw=0.5, ls='--')
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGURES / "flory_huggins_free_energy.png", dpi=200)
    plt.close(fig)
    print(f"  saved {FIGURES / 'flory_huggins_free_energy.png'}")


def plot_phase_diagram(N1=100, N2=1):
    """Phase diagram: χ vs φ with spinodal, binodal, and critical point."""
    phi_c, chi_c = critical_point(N1, N2)
    phi_sp, chi_sp = spinodal_curve(N1, N2)
    chi_bn, phi_bn_l, phi_bn_r = binodal_curve(N1, N2)

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.fill_betweenx(chi_sp, phi_sp[:len(chi_sp)],
                     phi_sp[::-1][:len(chi_sp)],
                     alpha=0.08, color='red', label='Spinodal region (unstable)')
    ax.plot(phi_sp, chi_sp, 'r--', lw=2, label='Spinodal')

    if len(chi_bn) > 0:
        ax.plot(phi_bn_l, chi_bn, 'b-', lw=2.5, label='Binodal')
        ax.plot(phi_bn_r, chi_bn, 'b-', lw=2.5)
        ax.fill_betweenx(chi_bn, phi_bn_l, phi_bn_r,
                         alpha=0.06, color='blue')

    ax.plot(phi_c, chi_c, 'ko', ms=10, zorder=5)
    ax.annotate(f'  Critical point\n  (φ={phi_c:.3f}, χ={chi_c:.3f})',
                xy=(phi_c, chi_c), fontsize=10,
                xytext=(phi_c + 0.15, chi_c + 0.02),
                arrowprops=dict(arrowstyle='->', color='black'))

    ax.set_xlabel("Volume fraction φ", fontsize=13)
    ax.set_ylabel("Flory–Huggins parameter χ  (∝ 1/T)", fontsize=13)
    ax.set_title(f"Phase diagram  (N₁ = {N1}, N₂ = {N2})", fontsize=14)
    ax.legend(fontsize=10, loc='upper right')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, chi_c * 3)
    ax.grid(True, alpha=0.3)

    ax2 = ax.twinx()
    chi_ticks = ax.get_yticks()
    chi_ticks = chi_ticks[chi_ticks > 0]
    ax2.set_ylim(ax.get_ylim())
    ax2.set_yticks(chi_ticks)
    ax2.set_yticklabels([f"{1/c:.1f}" for c in chi_ticks])
    ax2.set_ylabel("Effective temperature  T ∝ 1/χ  (a.u.)", fontsize=11)

    fig.tight_layout()
    fig.savefig(FIGURES / "flory_huggins_phase_diagram.png", dpi=200)
    plt.close(fig)
    print(f"  saved {FIGURES / 'flory_huggins_phase_diagram.png'}")


def plot_phase_diagram_temperature(N1=100, N2=1):
    """Phase diagram in T–φ coordinates (more physical)."""
    phi_c, chi_c = critical_point(N1, N2)
    T_c = 1.0 / chi_c

    phi_sp, chi_sp = spinodal_curve(N1, N2, n_points=300)
    T_sp = 1.0 / chi_sp

    chi_bn, phi_bn_l, phi_bn_r = binodal_curve(N1, N2, n_points=120)
    T_bn = 1.0 / chi_bn

    fig, ax = plt.subplots(figsize=(8, 6))

    if len(T_bn) > 0:
        ax.plot(phi_bn_l, T_bn, 'b-', lw=2.5, label='Binodal (coexistence)')
        ax.plot(phi_bn_r, T_bn, 'b-', lw=2.5)
        all_phi_bn = np.concatenate([phi_bn_l, phi_bn_r[::-1]])
        all_T_bn = np.concatenate([T_bn, T_bn[::-1]])
        ax.fill(all_phi_bn, all_T_bn, alpha=0.06, color='blue')

    # Reconstruct spinodal as left and right branches
    mid = len(phi_sp) // 2
    ax.plot(phi_sp[:mid], T_sp[:mid], 'r--', lw=2, label='Spinodal')
    ax.plot(phi_sp[mid:], T_sp[mid:], 'r--', lw=2)

    ax.plot(phi_c, T_c, 'ko', ms=10, zorder=5)
    ax.annotate(f'  UCST critical point\n  (φ={phi_c:.3f}, T*={T_c:.2f})',
                xy=(phi_c, T_c), fontsize=10,
                xytext=(phi_c + 0.12, T_c + T_c * 0.15),
                arrowprops=dict(arrowstyle='->', color='black'))

    ax.text(0.5, T_c * 0.5, 'Two-phase\nregion', fontsize=14,
            ha='center', va='center', color='navy', alpha=0.5)
    ax.text(0.5, T_c * 1.3, 'One phase\n(mixed)', fontsize=14,
            ha='center', va='center', color='green', alpha=0.5)

    ax.set_xlabel("Volume fraction φ", fontsize=13)
    ax.set_ylabel("Reduced temperature  T* = kT/ε", fontsize=13)
    ax.set_title(f"UCST phase diagram  (N₁ = {N1}, N₂ = {N2})", fontsize=14)
    ax.legend(fontsize=10, loc='upper right')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, T_c * 1.8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGURES / "flory_huggins_temperature_phase_diagram.png", dpi=200)
    plt.close(fig)
    print(f"  saved {FIGURES / 'flory_huggins_temperature_phase_diagram.png'}")


def plot_chain_length_effect():
    """Show how increasing chain length broadens the two-phase region."""
    fig, ax = plt.subplots(figsize=(8, 6))
    N_values = [1, 5, 20, 100, 1000]
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(N_values)))

    for N1, c in zip(N_values, colors):
        phi_c, chi_c = critical_point(N1, 1)
        phi_sp, chi_sp = spinodal_curve(N1, 1, n_points=300)
        T_sp = 1.0 / chi_sp
        T_c = 1.0 / chi_c

        mid = len(phi_sp) // 2
        ax.plot(phi_sp[:mid], T_sp[:mid] / T_c, '-', color=c, lw=2,
                label=f'N = {N1}')
        ax.plot(phi_sp[mid:], T_sp[mid:] / T_c, '-', color=c, lw=2)
        ax.plot(phi_c, 1.0, 'o', color=c, ms=8, zorder=5)

    ax.set_xlabel("Volume fraction φ", fontsize=13)
    ax.set_ylabel("T / T_c", fontsize=13)
    ax.set_title("Chain length widens the two-phase region", fontsize=14)
    ax.legend(fontsize=10, title="Polymer chain length")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.5)
    ax.axhline(1, color='grey', lw=0.5, ls=':')
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGURES / "flory_huggins_chain_length.png", dpi=200)
    plt.close(fig)
    print(f"  saved {FIGURES / 'flory_huggins_chain_length.png'}")


if __name__ == "__main__":
    print("Flory-Huggins phase separation")
    print("=" * 40)

    for N1 in [10, 100]:
        phi_c, chi_c = critical_point(N1, 1)
        print(f"  N₁={N1:>4d}: φ_c = {phi_c:.4f},  χ_c = {chi_c:.4f},  T_c* = {1/chi_c:.3f}")

    print("\nGenerating figures...")
    plot_free_energy()
    plot_phase_diagram()
    plot_phase_diagram_temperature()
    plot_chain_length_effect()
    print("Done.")
