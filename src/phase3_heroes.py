#!/usr/bin/env python3
r"""Phase 3 hero visualisations.

Stunning composite figures distilling the Phase 1 + Phase 2 results.

  * dispersion_surface_3d : the growth rate omega(k, phi0) as a 3D surface,
    with the spinodal contour (omega = 0 ridge) and the fastest-growing-mode
    crest highlighted.

Figures -> figures/.
"""

from __future__ import annotations

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path

from flory_huggins import d2f_dphi2
from nonlinear_analysis import set_style, ACCENT, INK

FIGURES = Path(__file__).resolve().parent.parent / "figures"
FIGURES.mkdir(exist_ok=True)


def dispersion_surface_3d(chi=1.5, N1=100, N2=1, kappa=0.5, M=1.0):
    """omega(k, phi0) as a 3D surface for the Flory-Huggins Cahn-Hilliard model."""
    set_style()
    phi = np.linspace(0.02, 0.66, 220)
    k = np.linspace(0.0, 2.2, 220)
    PHI, K = np.meshgrid(phi, k)
    FPP = d2f_dphi2(PHI, chi, N1, N2)
    W = -M * K**2 * (FPP + kappa * K**2)
    Wc = np.clip(W, -1.5, None)  # clip deep negatives for a clean surface

    fig = plt.figure(figsize=(11, 8))
    ax = fig.add_subplot(111, projection="3d")
    surf = ax.plot_surface(PHI, K, Wc, cmap="magma", rcount=120, ccount=120,
                           linewidth=0, antialiased=True, alpha=0.96,
                           vmin=-1.0, vmax=1.3)
    # zero plane (neutral stability)
    ax.plot_surface(PHI, K, np.zeros_like(Wc), color="grey", alpha=0.12,
                    linewidth=0)

    # fastest-growing-mode crest: k*(phi0) where f''<0
    kstar, wstar, pst = [], [], []
    for p in phi:
        fpp = d2f_dphi2(p, chi, N1, N2)
        if fpp < 0:
            ks = np.sqrt(-fpp / (2 * kappa))
            pst.append(p); kstar.append(ks)
            wstar.append(M * fpp**2 / (4 * kappa))
    ax.plot(pst, kstar, wstar, color="cyan", lw=3.2, zorder=10,
            label=r"fastest mode $k^*(\phi_0)$")

    ax.set_xlabel(r"$\phi_0$", labelpad=8)
    ax.set_ylabel(r"$k$", labelpad=8)
    ax.set_zlabel(r"$\omega(k,\phi_0)$", labelpad=6)
    ax.set_title("Dispersion surface of Cahn–Hilliard spinodal decomposition\n"
                 r"$\omega = -Mk^2\,(f''(\phi_0)+\kappa k^2)$  (Flory–Huggins, "
                 fr"$\chi={chi}$)", fontsize=13)
    ax.view_init(elev=26, azim=-58)
    ax.legend(loc="upper left", fontsize=10)
    fig.colorbar(surf, ax=ax, shrink=0.55, pad=0.02, label=r"$\omega$")
    fig.tight_layout()
    out = FIGURES / "phase3_dispersion_surface.png"
    fig.savefig(out, dpi=220)
    plt.close(fig)
    print(f"  saved {out}")


if __name__ == "__main__":
    print("Phase 3 hero visualisations")
    print("=" * 40)
    dispersion_surface_3d()
    print("Done.")
