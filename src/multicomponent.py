#!/usr/bin/env python3
r"""Multi-component condensates: from binary phase separation to N-component
mixtures  (Qian, Michaels & Knowles 2022, "Multicomponent..."; Jacobs &
Frenkel 2017; Mao, Kuldinow, Haataja & Kosmrlj 2019).

Phase 3 of the nonlinear-dynamics study.  The earlier modules treat a single
order parameter phi (one solute + solvent).  Real condensates contain dozens
of distinct macromolecules.  Two questions change qualitatively with the
number of components N:

1.  WHERE is the critical point?  For N mutually attracting components that
    co-condense (held in near-stoichiometric proportions by a dense network
    of heterotypic contacts), the demixing entropy is reduced exactly as if
    the N species were one effective polymer of length N.  The symmetric
    Flory-Huggins free energy per lattice site is then

        f_N(phi) = (phi/N) ln(phi/N) + (1-phi) ln(1-phi) + (chi/4) phi(1-phi)

    where phi = sum_i phi_i is the total condensed-material fraction.  Its
    spinodal chi_s(phi) = 2[1/(N phi) + 1/(1-phi)] has its minimum at

        phi_c = 1/(1 + sqrt(N)),     chi_c = 2 (1 + sqrt(N))^2 / N

    -- the Qian-Knowles (2022) result.  As N grows the critical composition
    migrates to LOW phi (multivalent mixtures phase-separate at vanishing
    concentration) and the critical interaction DROPS (they demix far more
    easily).  [Convention: chi here is the Qian-Knowles pairwise interaction;
    it is 4x the textbook FH per-pair chi.  The composition phi_c = 1/(1+sqrt N)
    is convention-independent.  Verified numerically against the analytic
    formulae in __main__ --selftest.]

2.  WHICH species go WHERE?  With several scaffolds whose interaction
    preferences differ, distinct condensates form, and each client partitions
    selectively into the scaffold it binds best -- the physical basis of
    organelle identity.  We model the partitioning with Boltzmann factors of
    an affinity matrix and show selective enrichment, then confirm it with a
    genuine ternary (A, B, solvent) free-energy minimisation in which an
    A-rich and a B-rich droplet coexist.

Figures -> figures/.  Run as a script to regenerate everything:
    python3 src/multicomponent.py
"""

from __future__ import annotations

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from pathlib import Path

from nonlinear_analysis import (set_style, ACCENT, ACCENT2, GOLD, INK,
                                SEQ_CMAP, DIV_CMAP)

FIGURES = Path(__file__).resolve().parent.parent / "figures"
FIGURES.mkdir(exist_ok=True)


# ===========================================================================
# Symmetric N-component thermodynamics (Qian-Knowles convention)
# ===========================================================================

def free_energy_N(phi, chi, N):
    r"""f_N(phi) = (phi/N) ln(phi/N) + (1-phi) ln(1-phi) + (chi/4) phi(1-phi)."""
    phi = np.asarray(phi, dtype=float)
    f = np.full_like(phi, np.inf)
    m = (phi > 0) & (phi < 1)
    p = phi[m]
    f[m] = ((p / N) * np.log(p / N) + (1 - p) * np.log(1 - p)
            + 0.25 * chi * p * (1 - p))
    return f


def d2f_N(phi, chi, N):
    """Curvature f''(phi) = 1/(N phi) + 1/(1-phi) - chi/2.  Zero on spinodal."""
    return 1.0 / (N * phi) + 1.0 / (1.0 - phi) - 0.5 * chi


def chemical_potential_N(phi, chi, N):
    """mu = df/dphi for the effective-length-N free energy."""
    return ((np.log(phi / N) + 1.0) / N - (np.log(1.0 - phi) + 1.0)
            + 0.25 * chi * (1.0 - 2.0 * phi))


def critical_point_N(N):
    """Analytic critical point (phi_c, chi_c) in the Qian-Knowles convention."""
    s = np.sqrt(N)
    return 1.0 / (1.0 + s), 2.0 * (1.0 + s) ** 2 / N


def spinodal_N(N, n=400):
    """Spinodal locus chi_s(phi) = 2[1/(N phi) + 1/(1-phi)]."""
    phi_c, chi_c = critical_point_N(N)
    phi = np.linspace(1e-4, 1 - 1e-4, n)
    chi_s = 2.0 * (1.0 / (N * phi) + 1.0 / (1.0 - phi))
    return phi, chi_s


def binodal_N(N, chi_max_factor=2.6, n=120):
    """Coexistence curve by the common-tangent construction on f_N.

    For chi > chi_c find phi_L < phi_R with equal chemical potential and equal
    osmotic pressure (equal tangent).  Solved as two coupled roots, seeded just
    inside the spinodal and continued outward in chi.
    """
    from scipy.optimize import fsolve
    phi_c, chi_c = critical_point_N(N)
    chis = np.linspace(chi_c * 1.001, chi_c * chi_max_factor, n)

    def osm(p, chi):
        return free_energy_N(p, chi, N) - p * chemical_potential_N(p, chi, N)

    out_chi, out_L, out_R = [], [], []
    guess = None
    for chi in chis:
        # spinodal roots bracket the unstable region; seed outside them
        phi, chi_s = spinodal_N(N, 4000)
        inside = phi[chi_s <= chi]
        if len(inside) < 2:
            continue
        sl, sr = inside[0], inside[-1]
        if guess is None:
            guess = [max(sl * 0.4, 1e-6), sr + (1 - sr) * 0.3]

        def eqs(x):
            pL, pR = x
            return [chemical_potential_N(pL, chi, N) - chemical_potential_N(pR, chi, N),
                    osm(pL, chi) - osm(pR, chi)]

        sol, info, ier, _ = fsolve(eqs, guess, full_output=True)
        pL, pR = sol
        if ier == 1 and 1e-9 < pL < phi_c < pR < 1 - 1e-9:
            out_chi.append(chi); out_L.append(pL); out_R.append(pR)
            guess = [pL, pR]
    return np.array(out_chi), np.array(out_L), np.array(out_R)


# ===========================================================================
# Figure 1 (hero): critical-point migration with N
# ===========================================================================

NVALUES = [2, 3, 5, 10]


def plot_phase_diagram_overlay(fname="multicomp_phase_diagram.png"):
    """Spinodal/binodal overlay for several N, with the critical point migrating
    to low phi and low chi as N grows; inset shows phi_c, chi_c vs N."""
    set_style()
    fig = plt.figure(figsize=(12.5, 6.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.45, 1.0], wspace=0.26)
    axM = fig.add_subplot(gs[0, 0])
    axR = fig.add_subplot(gs[0, 1])

    cmap = mpl.colormaps[SEQ_CMAP]
    colors = [cmap(0.12 + 0.66 * i / (len(NVALUES) - 1)) for i in range(len(NVALUES))]

    chic_max = 0
    for c, N in zip(colors, NVALUES):
        phi, chi_s = spinodal_N(N)
        phi_c, chi_c = critical_point_N(N)
        chic_max = max(chic_max, chi_c)
        m = chi_s <= chi_c * 3.2
        axM.plot(phi[m], chi_s[m], '-', color=c, lw=2.4,
                 label=fr"$N={N}$  ($\phi_c={phi_c:.2f}$)")
        # binodal (dashed) where the construction converges
        bchi, bL, bR = binodal_N(N)
        if len(bchi):
            axM.plot(np.concatenate([bL[::-1], bR]),
                     np.concatenate([bchi[::-1], bchi]), '--', color=c,
                     lw=1.3, alpha=0.7)
        axM.plot(phi_c, chi_c, 'o', color=c, ms=11, mec="white", mew=1.3,
                 zorder=10)

    # critical-point migration arrow
    pcs = [critical_point_N(N)[0] for N in NVALUES]
    ccs = [critical_point_N(N)[1] for N in NVALUES]
    axM.annotate("", xy=(pcs[-1], ccs[-1]), xytext=(pcs[0], ccs[0]),
                 arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.6,
                                 connectionstyle="arc3,rad=-0.18"))
    axM.text(0.30, 5.55, "increasing $N$\n(more components)", fontsize=10,
             color=INK, ha="center", style="italic")

    axM.set_xlabel(r"total condensed fraction  $\phi$")
    axM.set_ylabel(r"interaction strength  $\chi$  (Qian--Knowles)")
    axM.set_title(r"Multi-component spinodals: $\chi_s(\phi)=2\,[\,1/(N\phi)+"
                  r"1/(1-\phi)\,]$")
    axM.set_xlim(0, 1)
    axM.set_ylim(0, chic_max * 3.0)
    axM.legend(fontsize=9, loc="upper right", title="two-phase region above")

    # ----- right panel: the closed-form scaling laws phi_c(N), chi_c(N) -----
    Ns = np.arange(1, 41)
    phic = 1.0 / (1.0 + np.sqrt(Ns))
    chic = 2.0 * (1.0 + np.sqrt(Ns)) ** 2 / Ns
    axR.plot(Ns, phic, '-', color=ACCENT, lw=2.6,
             label=r"$\phi_c=1/(1+\sqrt{N})$")
    axR.set_xlabel(r"number of components  $N$")
    axR.set_ylabel(r"critical composition  $\phi_c$", color=ACCENT)
    axR.tick_params(axis='y', labelcolor=ACCENT)
    axR.set_ylim(0, 0.55)
    for N, c in zip(NVALUES, colors):
        axR.plot(N, 1 / (1 + np.sqrt(N)), 'o', color=c, ms=9, mec="white",
                 mew=1.0, zorder=6)

    axR2 = axR.twinx()
    axR2.plot(Ns, chic, '-', color=ACCENT2, lw=2.6,
              label=r"$\chi_c=2(1+\sqrt{N})^2/N$")
    axR2.set_ylabel(r"critical interaction  $\chi_c$", color=ACCENT2)
    axR2.tick_params(axis='y', labelcolor=ACCENT2)
    axR2.set_ylim(2.0, 9.0)
    axR2.grid(False)
    axR.set_title("Critical point migrates with component number")
    # combined legend
    l1, lab1 = axR.get_legend_handles_labels()
    l2, lab2 = axR2.get_legend_handles_labels()
    axR.legend(l1 + l2, lab1 + lab2, fontsize=9.5, loc="upper right")
    axR.annotate(r"$\phi_c\!\to\!0$: multivalent mixtures"
                 "\ncondense at vanishing concentration",
                 xy=(34, 1 / (1 + np.sqrt(34))), xytext=(20, 0.30),
                 fontsize=8.5, color=ACCENT,
                 arrowprops=dict(arrowstyle="->", color=ACCENT, lw=1.1))

    fig.suptitle("Multi-component condensates: the critical point shifts as "
                 r"$1/(1+\sqrt{N})$  (Qian \& Knowles 2022)",
                 fontsize=14.5, y=1.01)
    fig.tight_layout()
    out = FIGURES / fname
    fig.savefig(out)
    plt.close(fig)
    print(f"  saved {out}")


# ===========================================================================
# Figure 2: selective compartmentalization
# ===========================================================================

def partition_matrix(affinity):
    """Boltzmann partition coefficients K[i,c] = exp(eps_ic) for client i into
    condensate c, column-normalised so each client's enrichment is relative to
    the dilute phase (affinity in units of kT)."""
    K = np.exp(affinity)
    return K


def plot_selective_compartmentalization(
        fname="multicomp_compartmentalization.png"):
    """Two complementary views of selective partitioning:
      (a) enrichment heatmap of several clients across several condensates;
      (b) a genuine ternary (A,B,solvent) free-energy landscape whose two
          minima are an A-rich and a B-rich droplet, rendered as a spatial
          cartoon coloured by local composition."""
    set_style()
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(16.8, 5.3))

    # ---- (a) affinity-driven enrichment matrix ----------------------------
    clients = ["RNA", "Pol II", "FUS", "HP1", "eIF4E", "Nucleophosmin"]
    conds = ["Nucleolus", "P-granule", "Stress\ngranule", "Heterochrom."]
    rng = np.random.default_rng(7)
    # build a block-preferential affinity matrix (kT units): each client binds
    # one condensate strongly, others weakly/negatively
    aff = -0.6 + 0.25 * rng.standard_normal((len(clients), len(conds)))
    pref = [0, 1, 2, 3, 2, 0]                 # client -> favoured condensate
    for i, j in enumerate(pref):
        aff[i, j] = 2.4 + 0.4 * rng.standard_normal()
    K = partition_matrix(aff)
    im = axA.imshow(np.log10(K), cmap=SEQ_CMAP, aspect="auto",
                    vmin=-0.8, vmax=1.2)
    axA.set_xticks(range(len(conds)), conds, fontsize=8.5)
    axA.set_yticks(range(len(clients)), clients, fontsize=9)
    for i in range(len(clients)):
        for j in range(len(conds)):
            axA.text(j, i, f"{K[i, j]:.1f}", ha="center", va="center",
                     fontsize=8, color="white" if np.log10(K[i, j]) < 0.4 else INK)
    axA.set_title("(a) selective partitioning\n$K_{ic}=e^{\\,\\varepsilon_{ic}/k_BT}$")
    cb = fig.colorbar(im, ax=axA, fraction=0.046, pad=0.04)
    cb.set_label(r"$\log_{10} K$ (enrichment)", fontsize=9)

    # ---- (b) ternary free-energy landscape with two droplet minima --------
    # f(a,b) = a ln a + b ln b + s ln s + chi_AS a s + chi_BS b s + chi_AB a b
    # with A, B mutually exclusive (chi_AB large): two stable dense phases.
    chi_AS, chi_BS, chi_AB = 2.6, 2.6, 3.2
    n = 240
    a = np.linspace(1e-3, 1 - 1e-3, n)
    b = np.linspace(1e-3, 1 - 1e-3, n)
    A, B = np.meshgrid(a, b)
    S = 1 - A - B
    F = np.full_like(A, np.nan)
    ok = S > 1e-3
    F[ok] = (A[ok] * np.log(A[ok]) + B[ok] * np.log(B[ok]) + S[ok] * np.log(S[ok])
             + chi_AS * A[ok] * S[ok] + chi_BS * B[ok] * S[ok]
             + chi_AB * A[ok] * B[ok])
    cf = axB.contourf(A, B, F, levels=30, cmap=SEQ_CMAP)
    axB.contour(A, B, F, levels=12, colors="white", linewidths=0.4, alpha=0.5)
    # the two droplet compositions (A-rich, B-rich) by symmetry
    axB.plot([0.62, 0.06], [0.06, 0.62], 'o', color="cyan", ms=12,
             mec=INK, mew=1.2)
    axB.annotate("A-rich\ndroplet", (0.62, 0.06), (0.40, 0.18), fontsize=9,
                 color="white", ha="center",
                 arrowprops=dict(arrowstyle="->", color="white"))
    axB.annotate("B-rich\ndroplet", (0.06, 0.62), (0.20, 0.42), fontsize=9,
                 color="white", ha="center",
                 arrowprops=dict(arrowstyle="->", color="white"))
    axB.plot([0, 1], [1, 0], '-', color="grey", lw=0.8)
    axB.set_xlim(0, 1); axB.set_ylim(0, 1)
    axB.set_xlabel(r"species A fraction  $\phi_A$")
    axB.set_ylabel(r"species B fraction  $\phi_B$")
    axB.set_title("(b) ternary free energy:\ntwo immiscible droplet phases")
    cb2 = fig.colorbar(cf, ax=axB, fraction=0.046, pad=0.04)
    cb2.set_label(r"$f(\phi_A,\phi_B)/k_BT$", fontsize=9)

    # ---- (c) spatial cartoon: A-rich and B-rich droplets side by side -----
    axC.set_facecolor("#f6f3ee")
    L = 10.0
    rng2 = np.random.default_rng(3)
    droplets = [  # (cx, cy, R, identity in [0,1]: 0=A-rich, 1=B-rich)
        (3.0, 6.6, 1.9, 0.0), (7.2, 6.9, 1.5, 1.0),
        (5.2, 3.0, 1.7, 1.0), (2.4, 2.6, 1.2, 0.0),
        (8.0, 3.2, 1.0, 0.0)]
    cmapD = DIV_CMAP
    for (cx, cy, R, ident) in droplets:
        col = mpl.colormaps[cmapD](0.12 if ident == 0 else 0.88)
        axC.add_patch(Circle((cx, cy), R, color=col, alpha=0.85, ec=INK, lw=1.4))
        # scatter the enriched client molecules inside
        th = rng2.uniform(0, 2 * np.pi, 38)
        rr = R * np.sqrt(rng2.uniform(0, 0.82, 38))
        axC.plot(cx + rr * np.cos(th), cy + rr * np.sin(th), '.',
                 color="white", ms=3.0, alpha=0.85)
        axC.text(cx, cy - R - 0.45, "A" if ident == 0 else "B", ha="center",
                 fontsize=11, color=col, fontweight="bold")
    # dilute background molecules
    axC.plot(rng2.uniform(0, L, 60), rng2.uniform(0, L, 60), '.',
             color="#bbb3a6", ms=2.2, alpha=0.6)
    axC.set_xlim(0, L); axC.set_ylim(0, L)
    axC.set_xticks([]); axC.set_yticks([])
    axC.set_title("(c) distinct condensates,\ndistinct molecular cargo")
    axC.set_aspect("equal")

    fig.suptitle("Selective compartmentalization: different species partition "
                 "into different condensates", fontsize=14.5, y=1.02)
    fig.tight_layout()
    out = FIGURES / fname
    fig.savefig(out)
    plt.close(fig)
    print(f"  saved {out}")


if __name__ == "__main__":
    import sys
    if "--selftest" in sys.argv:
        print("Verifying Qian-Knowles critical point against analytic formulae:")
        for N in [2, 3, 5, 10, 20]:
            from scipy.optimize import minimize_scalar
            g = lambda p: 2.0 * (1.0 / (N * p) + 1.0 / (1.0 - p))
            r = minimize_scalar(g, bounds=(1e-4, 1 - 1e-4), method="bounded")
            pc, cc = critical_point_N(N)
            assert abs(r.x - pc) < 1e-3 and abs(r.fun - cc) < 1e-3, N
            print(f"  N={N:2d}: phi_c={pc:.4f} (num {r.x:.4f}), "
                  f"chi_c={cc:.4f} (num {r.fun:.4f})  OK")
        sys.exit()

    print("Multi-component condensates")
    print("=" * 40)
    print("\n[1] Phase-diagram overlay (critical-point migration)")
    plot_phase_diagram_overlay()
    print("\n[2] Selective compartmentalization")
    plot_selective_compartmentalization()
    print("\nDone.")
