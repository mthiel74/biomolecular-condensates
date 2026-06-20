#!/usr/bin/env python3
r"""Flory-Huggins extensions: the stickers-and-spacers model.

Phase 2 thermodynamics.  Multivalent biomolecules phase-separate far more
readily than simple FH theory predicts, because their interactions are
specific and saturable ("stickers" joined by flexible "spacers";
Harmon & Pappu 2017).  We add a saturable association free energy to the
symmetric Flory-Huggins free energy and show how increasing sticker valence
expands the two-phase region and lowers the critical interaction.

Association free energy (Semenov-Rubinstein sticky-polymer theory): with
sticker volume fraction phi_s = v*phi (valence v) and association constant
K_a, the equilibrium bond fraction p solves the law of mass action

    p / (1-p)^2 = K_a phi_s        =>    p = [(2x+1) - sqrt(4x+1)] / (2x),
                                          x = K_a phi_s,

and the bonding free energy density is

    f_bond(phi) = phi_s [ ln(1-p) + p/2 ].

This is bounded (unlike a naive -rho ln(1+K rho)), so it produces a modest,
physical shift rather than diverging.  Adding f_bond to the symmetric FH
free energy lowers the critical chi and shifts phi_c -- multivalency drives
liquid-liquid phase separation.

Figure -> figures/fh_stickers_phase_diagram.png
"""

from __future__ import annotations

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path

from nonlinear_analysis import set_style, ACCENT, ACCENT2, GOLD, INK

FIGURES = Path(__file__).resolve().parent.parent / "figures"
FIGURES.mkdir(exist_ok=True)

N1 = N2 = 1  # symmetric (protein-like) baseline; FH critical chi = 2


def fh_symmetric(phi, chi):
    p = np.clip(phi, 1e-12, 1 - 1e-12)
    return p * np.log(p) / N1 + (1 - p) * np.log(1 - p) / N2 + chi * p * (1 - p)


def bond_fraction(phi_s, Ka):
    x = Ka * np.maximum(phi_s, 1e-12)
    return ((2 * x + 1) - np.sqrt(4 * x + 1)) / (2 * np.maximum(x, 1e-12))


def f_bond(phi, v, Ka):
    if Ka <= 0 or v <= 0:
        return np.zeros_like(np.asarray(phi, float))
    phi_s = v * phi
    p = bond_fraction(phi_s, Ka)
    return phi_s * (np.log(np.clip(1 - p, 1e-12, 1)) + p / 2)


def f_total(phi, chi, v, Ka):
    return fh_symmetric(phi, chi) + f_bond(phi, v, Ka)


def _fbond_pp(phi, v, Ka, h=1e-4):
    return (f_bond(phi + h, v, Ka) - 2 * f_bond(phi, v, Ka)
            + f_bond(phi - h, v, Ka)) / h**2


def spinodal_sticker(v, Ka, n=4000):
    """chi on the spinodal as a function of phi (f''_total = 0)."""
    phi = np.linspace(0.005, 0.995, n)
    extra = _fbond_pp(phi, v, Ka) if Ka > 0 else 0.0
    chi = 0.5 * (1 / (N1 * phi) + 1 / (N2 * (1 - phi)) + extra)
    return phi, chi


def critical_point_sticker(v, Ka):
    phi, chi = spinodal_sticker(v, Ka, n=8000)
    i = int(np.argmin(chi))
    return phi[i], chi[i]


def binodal_hull(chi, v, Ka, ngrid=4000):
    """Binodal (coexistence) endpoints from the lower convex hull of f_total."""
    phi = np.linspace(1e-4, 1 - 1e-4, ngrid)
    f = f_total(phi, chi, v, Ka)
    ok = np.isfinite(f)
    phi, f = phi[ok], f[ok]
    hull = []
    for x, y in zip(phi, f):
        while len(hull) >= 2:
            (x1, y1), (x2, y2) = hull[-2], hull[-1]
            if (x2 - x1) * (y - y1) - (y2 - y1) * (x - x1) <= 0:
                hull.pop()
            else:
                break
        hull.append((x, y))
    hx = np.array([h[0] for h in hull])
    j = int(np.argmax(np.diff(hx)))
    return hx[j], hx[j + 1]


def binodal_curve_sticker(v, Ka, n=60):
    _, chi_c = critical_point_sticker(v, Ka)
    chis = np.linspace(chi_c + 1e-3, chi_c + 2.2, n)
    L, R, C = [], [], []
    for chi in chis:
        pl, pr = binodal_hull(chi, v, Ka)
        if pr - pl > 1e-3:
            L.append(pl); R.append(pr); C.append(chi)
    return np.array(C), np.array(L), np.array(R)


def plot_sticker_phase_diagram():
    set_style()
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(16.5, 5.2))
    cases = [(1, 0.0, "no stickers (FH)"), (2, 2.0, "v=2, $K_a$=2"),
             (3, 2.0, "v=3, $K_a$=2"), (4, 4.0, "v=4, $K_a$=4")]
    cmap = mpl.colormaps["plasma"]
    colors = [cmap(0.1 + 0.75 * i / (len(cases) - 1)) for i in range(len(cases))]

    # (a) free energy at a fixed chi
    chi_demo = 1.6
    phi = np.linspace(1e-3, 1 - 1e-3, 500)
    for (v, Ka, lab), c in zip(cases, colors):
        axA.plot(phi, f_total(phi, chi_demo, v, Ka), color=c, lw=2.2, label=lab)
    axA.set_xlabel(r"volume fraction  $\phi$")
    axA.set_ylabel(r"free energy  $f(\phi)$")
    axA.set_title(fr"(a) free energy at $\chi={chi_demo}$")
    axA.legend(fontsize=9)

    # (b) spinodal curves + critical points
    for (v, Ka, lab), c in zip(cases, colors):
        phi_s, chi_s = spinodal_sticker(v, Ka)
        m = chi_s < 4.5
        axB.plot(phi_s[m], chi_s[m], color=c, lw=2.2, label=lab)
        pc, cc = critical_point_sticker(v, Ka)
        axB.plot(pc, cc, 'o', color=c, ms=8, mec="white", mew=0.6, zorder=5)
    axB.set_xlabel(r"volume fraction  $\phi$")
    axB.set_ylabel(r"interaction  $\chi$")
    axB.set_title("(b) stickers expand the two-phase region")
    axB.set_ylim(0, 3.2)
    axB.legend(fontsize=9)
    axB.text(0.5, 0.4, "two-phase below each curve", ha="center", fontsize=9,
             color=INK, alpha=0.7)

    # (c) critical chi vs valence (the phase-separation propensity)
    valences = np.arange(1, 9)
    for Ka, c in [(1.0, ACCENT2), (2.0, ACCENT), (4.0, GOLD)]:
        ccs = [critical_point_sticker(v if v > 1 else 1, Ka if v > 1 else 0)[1]
               for v in valences]
        axC.plot(valences, ccs, 'o-', color=c, lw=2, ms=6, label=fr"$K_a={Ka:g}$")
    axC.axhline(2.0, color=INK, ls=":", lw=1, label="FH baseline ($\\chi_c=2$)")
    axC.set_xlabel(r"sticker valence  $v$")
    axC.set_ylabel(r"critical interaction  $\chi_c$")
    axC.set_title("(c) multivalency lowers the LLPS threshold")
    axC.legend(fontsize=9)

    fig.suptitle("Stickers-and-spacers: how multivalency drives liquid-liquid "
                 "phase separation (Harmon & Pappu 2017)", fontsize=14.5, y=1.01)
    fig.tight_layout()
    out = FIGURES / "fh_stickers_phase_diagram.png"
    fig.savefig(out); plt.close(fig)
    print(f"  saved {out}")


if __name__ == "__main__":
    print("Flory-Huggins stickers-and-spacers extension")
    print("=" * 44)
    for v, Ka in [(1, 0), (2, 2), (4, 4)]:
        pc, cc = critical_point_sticker(v, Ka)
        print(f"  v={v}, Ka={Ka}: phi_c={pc:.3f}, chi_c={cc:.3f}")
    print("\nGenerating figure...")
    plot_sticker_phase_diagram()
    print("Done.")