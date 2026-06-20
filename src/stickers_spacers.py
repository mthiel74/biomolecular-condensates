#!/usr/bin/env python3
r"""Stickers-and-spacers: how molecular valence shapes the phase diagram
(Harmon, Holehouse, Rosen & Pappu 2017, eLife; Choi, Holehouse & Pappu 2020,
Annu. Rev. Biophys.).

Phase 2 of the nonlinear-dynamics study.  The abstract Flory interaction
parameter chi of Phase 1 hides the molecular architecture of the
condensate-forming protein.  The stickers-and-spacers picture makes it
concrete: multivalent proteins phase separate through *specific, saturable*
associative contacts ("stickers") connected by flexible "spacers".  We add a
reversible-association free energy to the Flory-Huggins functional:

    f(phi) = phi ln(phi)/N1 + (1-phi) ln(1-phi)/N2 + chi phi(1-phi)
             + f_assoc(phi)
    f_assoc(phi) = - rho_s ln(1 + K_assoc rho_s),   rho_s = (valence/N1) phi

rho_s is the density of stickers (valence per chain times the chain density
phi/N1), K_assoc is the sticker-sticker association constant.  The term is
negative and grows with density, so association lowers the free energy of the
dense phase and *promotes* condensation: increasing the valence or the
association strength widens the two-phase region and shifts the binodal to
lower chi (i.e. weaker mean-field attraction is needed to demix).

Binodals are computed by the parameter-free lower-convex-hull (double-tangent)
construction, which is robust for any free energy -- in particular it does not
suffer the fsolve failure that the 2-D common-tangent solver hits for the
strongly asymmetric N1=100 case.

Figures -> figures/.  Run as a script:  python3 src/stickers_spacers.py
"""

from __future__ import annotations

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path

from flory_huggins import free_energy_density, critical_point
from nonlinear_analysis import set_style, ACCENT, ACCENT2, GOLD, INK

FIGURES = Path(__file__).resolve().parent.parent / "figures"
FIGURES.mkdir(exist_ok=True)


# ===========================================================================
# Free energy with reversible association
# ===========================================================================

def association_free_energy(phi, valence, K_assoc, N1=100):
    r"""Reversible-association contribution  -rho_s ln(1 + K rho_s).

    rho_s = (valence / N1) * phi  is the sticker density (stickers per chain
    times chain volume fraction).  valence = 0 recovers plain Flory-Huggins.
    """
    rho_s = (valence / N1) * np.asarray(phi, dtype=float)
    return -rho_s * np.log1p(K_assoc * rho_s)


def free_energy_sticker(phi, chi, valence=0, K_assoc=0.0, N1=100, N2=1):
    """Total Flory-Huggins + association free energy density."""
    f = free_energy_density(phi, chi, N1, N2)
    if valence > 0 and K_assoc > 0:
        f = f + association_free_energy(phi, valence, K_assoc, N1)
    return f


# ===========================================================================
# Robust binodal via the lower convex hull (double-tangent construction)
# ===========================================================================

def binodal_endpoints(f_vals, phi, atol=1e-4):
    r"""Endpoints of the two-phase region from the lower convex envelope of f.

    Where the convex hull of (phi, f) lies strictly below the curve f, the
    system lowers its free energy by demixing along a tie-line; the segment's
    endpoints are the coexisting compositions (phi_L, phi_R).  Returns the
    widest such tie-line (the outer binodal), or None if f is convex.
    """
    n = len(phi)
    # lower convex hull of the point set (monotone chain on sorted phi)
    hull = []
    for i in range(n):
        while len(hull) >= 2:
            x1, y1 = phi[hull[-2]], f_vals[hull[-2]]
            x2, y2 = phi[hull[-1]], f_vals[hull[-1]]
            x3, y3 = phi[i], f_vals[i]
            # pop if last point is above the line (x1->x3): not on lower hull
            if (y2 - y1) * (x3 - x1) >= (y3 - y1) * (x2 - x1):
                hull.pop()
            else:
                break
        hull.append(i)
    # find the longest hull edge that skips interior points (= tie-line)
    best = None
    best_gap = 0.0
    for a, b in zip(hull[:-1], hull[1:]):
        if b - a > 1 and (phi[b] - phi[a]) > best_gap:
            # confirm the curve bulges above this chord somewhere between
            xs = phi[a:b + 1]
            chord = f_vals[a] + (f_vals[b] - f_vals[a]) * \
                (xs - phi[a]) / (phi[b] - phi[a])
            if np.max(f_vals[a:b + 1] - chord) > atol:
                best_gap = phi[b] - phi[a]
                best = (phi[a], phi[b])
    return best


def sticker_binodal(valence, K_assoc, N1=100, N2=1, n_chi=140, n_phi=4000):
    """Trace the binodal phi_L(chi), phi_R(chi) for a sticker free energy."""
    phi = np.linspace(1e-5, 1 - 1e-5, n_phi)
    phi_c0, chi_c0 = critical_point(N1, N2)
    chi_vals = np.linspace(0.2 * chi_c0, 3.0 * chi_c0, n_chi)
    out_chi, out_l, out_r = [], [], []
    for chi in chi_vals:
        f = free_energy_sticker(phi, chi, valence, K_assoc, N1, N2)
        ends = binodal_endpoints(f, phi)
        if ends is not None:
            out_chi.append(chi)
            out_l.append(ends[0])
            out_r.append(ends[1])
    return np.array(out_chi), np.array(out_l), np.array(out_r)


# ===========================================================================
# Figures
# ===========================================================================

def plot_sticker_phase_diagram(N1=100, N2=1, K_assoc=80.0,
                               valences=(0, 2, 4, 8),
                               fname="sticker_phase_diagram.png"):
    """Overlay binodals: plain FH vs sticker-modified for increasing valence."""
    set_style()
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.5, 5.6))
    cmap = mpl.colormaps["plasma"]

    # (a) binodal overlay
    for j, v in enumerate(valences):
        chi, pL, pR = sticker_binodal(v, K_assoc, N1, N2)
        if len(chi) == 0:
            continue
        c = INK if v == 0 else cmap(0.15 + 0.7 * j / max(1, len(valences) - 1))
        lab = "plain FH (valence 0)" if v == 0 else fr"valence $= {v}$"
        axA.plot(pL, chi, color=c, lw=2.4, label=lab)
        axA.plot(pR, chi, color=c, lw=2.4)
        # shade the two-phase region of the highest valence
        if v == max(valences):
            axA.fill_betweenx(chi, pL, pR, color=c, alpha=0.08)
    phi_c0, chi_c0 = critical_point(N1, N2)
    axA.axhline(chi_c0, color="grey", lw=0.8, ls=":")
    axA.text(0.62, chi_c0 * 1.03, r"plain-FH $\chi_c$", color="grey",
             fontsize=9)
    axA.set_xlabel(r"volume fraction  $\phi$")
    axA.set_ylabel(r"Flory parameter  $\chi \propto 1/T$")
    axA.set_title("(a) valence widens the two-phase region\n"
                  "and lowers the demixing $\\chi$")
    axA.legend(fontsize=9, loc="upper right")
    axA.set_xlim(0, 1)
    axA.set_ylim(0, chi_c0 * 3)
    axA.annotate("", xy=(0.5, chi_c0 * 0.55), xytext=(0.5, chi_c0 * 1.6),
                 arrowprops=dict(arrowstyle="->", color=ACCENT, lw=1.6))
    axA.text(0.515, chi_c0 * 0.8, "stickers demix\nat higher $T$\n(lower $\\chi$)",
             color=ACCENT, fontsize=9, va="center")

    # (b) the saturating association free energy that drives it
    phi = np.linspace(1e-4, 1 - 1e-4, 500)
    for j, v in enumerate(valences):
        if v == 0:
            continue
        c = cmap(0.15 + 0.7 * j / max(1, len(valences) - 1))
        fa = association_free_energy(phi, v, K_assoc, N1)
        axB.plot(phi, fa, color=c, lw=2.2, label=fr"valence $= {v}$")
    axB.set_xlabel(r"volume fraction  $\phi$")
    axB.set_ylabel(r"association free energy  $f_{\rm assoc}(\phi)$")
    axB.set_title(r"(b) $f_{\rm assoc}=-\rho_s\ln(1+K\rho_s)$,"
                  "\nsaturable & valence-scaled")
    axB.legend(fontsize=9, loc="lower left")
    axB.set_xlim(0, 1)

    fig.suptitle("Stickers-and-spacers: molecular valence reshapes the "
                 "condensate phase diagram  (Pappu group)",
                 fontsize=14.5, y=1.02)
    fig.tight_layout()
    out = FIGURES / fname
    fig.savefig(out)
    plt.close(fig)
    print(f"  saved {out}")


def report_widening(N1=100, K_assoc=80.0, valences=(0, 2, 4, 8)):
    """Print the demixing-onset chi (lowest chi with a binodal) per valence."""
    print("  valence   chi_onset   phi_dense(max-chi)")
    for v in valences:
        chi, pL, pR = sticker_binodal(v, K_assoc, N1)
        if len(chi):
            print(f"  {v:>5}     {chi.min():.4f}      {pR.max():.3f}")
        else:
            print(f"  {v:>5}     (no demixing)")


if __name__ == "__main__":
    print("Stickers-and-spacers phase behaviour")
    print("=" * 40)
    # validate the convex-hull binodal against the known FH critical point
    phi_c, chi_c = critical_point(100, 1)
    print(f"  plain-FH critical point: phi_c={phi_c:.4f}, chi_c={chi_c:.4f}")
    report_widening()
    print("\nGenerating figure...")
    plot_sticker_phase_diagram()
    print("Done.")
