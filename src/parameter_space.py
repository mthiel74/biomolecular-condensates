#!/usr/bin/env python3
r"""Full parameter-space exploration: a map of every regime the condensate
model can access.

Phase 3 of the nonlinear-dynamics study; the synthesis.  The preceding modules
each lived at one point of a much larger parameter space.  Here we sweep four
2D control planes and draw the regime boundaries -- spinodal vs nucleation vs
homogeneous, dissolution vs emulsion vs macrophase, one-phase vs two-phase,
and the rate-enhancement landscape -- locating the tricritical points and
codimension-2 corners where three regimes meet.

Planes
------
A.  (phi_0, kappa)    : thermodynamic stability of the double well, with the
                        gradient cost kappa setting a finite-size cutoff.
                        Spinodal |phi_0| < sqrt(a/3b); binodal |phi_0| <
                        sqrt(a/b); inside the spinodal a pattern of wavelength
                        lambda* = 2 pi sqrt(2 kappa/|f''|) grows -- if it fits
                        the box.  Regimes: homogeneous / nucleation / spinodal.

B.  (k_fuel, lambda)  : the ACTIVE plane (Active Model B+ with a chemical
                        reaction).  Short 2D spectral runs classify each point
                        as dissolution / microphase (arrested emulsion) /
                        macrophase by the order parameter and the domain length.
                        A codimension-2 corner joins all three.

C.  (chi, N)          : Flory-Huggins.  Two-phase iff chi > chi_c(N) =
                        (1/2)(1 + 1/sqrt N)^2; the critical line chi_c(N) is a
                        line of critical points migrating to low chi with N.

D.  (K, p)            : mass-conserving rate enhancement k_eff/k_out over the
                        partition coefficient K and dense area fraction p, with
                        the optimal ridge p* ~ 1/K (the honest interior maximum,
                        not the naive K^2).

Grand summary figure -> figures/.  Run as a script:
    python3 src/parameter_space.py
    python3 src/parameter_space.py --quick
"""

from __future__ import annotations

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, BoundaryNorm, ListedColormap
from pathlib import Path

from nonlinear_analysis import (set_style, ACCENT, ACCENT2, GOLD, INK,
                                effective_rate_multiplier)
from flory_huggins import critical_point
from active_models import ActiveCH2D

FIGURES = Path(__file__).resolve().parent.parent / "figures"
FIGURES.mkdir(exist_ok=True)

# regime palette (shared)
REG_COLORS = ["#e9e4d8", "#7fb3d5", "#c81e5b", "#e8a33d"]   # homog/micro/macro/...


# ===========================================================================
# Plane A:  (phi_0, kappa)  -- spinodal / nucleation / homogeneous (analytic)
# ===========================================================================

def plane_phi0_kappa(a=1.0, b=1.0, M=1.0, box=20.0, n=320):
    """Regime field over (phi_0, kappa) for the symmetric double well.

    0 homogeneous (stable, outside binodal, or pattern too large for the box)
    1 nucleation  (metastable: between binodal and spinodal)
    2 spinodal    (unstable AND lambda* < box -> spontaneous decomposition)
    """
    phi0 = np.linspace(-1.0, 1.0, n)
    kappa = np.linspace(0.05, 8.0, n)
    PHI, KAP = np.meshgrid(phi0, kappa)
    fpp = -a + 3 * b * PHI**2                       # f''(phi0), double well
    phi_s = np.sqrt(a / (3 * b))                    # spinodal edge
    phi_b = np.sqrt(a / b)                          # binodal edge (|phi|<phi_b)
    reg = np.zeros_like(PHI)
    metastable = (np.abs(PHI) < phi_b) & (fpp >= 0)
    reg[metastable] = 1
    unstable = fpp < 0
    lam_star = np.full_like(PHI, np.inf)
    lam_star[unstable] = 2 * np.pi * np.sqrt(2 * KAP[unstable] / (-fpp[unstable]))
    reg[unstable & (lam_star < box)] = 2
    reg[unstable & (lam_star >= box)] = 1           # unstable but pattern > box
    return dict(phi0=phi0, kappa=kappa, reg=reg, phi_s=phi_s, phi_b=phi_b,
                lam_star=lam_star, box=box)


# ===========================================================================
# Plane B:  (k_fuel, lambda)  -- active regimes (simulation sweep)
# ===========================================================================

def plane_fuel_activity(nk=11, nl=9, N=64, L=64.0, T=260.0, dt=0.05,
                        kmax=0.30, lmax=2.4, seed=4, verbose=True):
    """Short 2D Active-Model-B+ runs over (k_react, lambda); classify each by
    the order parameter std(phi) and the structure-factor domain length.

    0 dissolution (std small -> homogeneous)
    1 microphase  (ordered but short length scale -> arrested emulsion)
    2 macrophase  (ordered, long length scale -> coarsened domains)
    """
    ks = np.linspace(0.0, kmax, nk)
    lams = np.linspace(0.0, lmax, nl)
    op = np.zeros((nl, nk))
    Ld = np.zeros((nl, nk))
    n_steps = int(T / dt)
    for i, lam in enumerate(lams):
        row = []
        for j, kr in enumerate(ks):
            sim = ActiveCH2D(L=L, N=N, dt=dt, k_react=kr, phi_ss=0.0, lam=lam,
                             A=2.6)
            sim.initialize(phi_mean=0.0, noise_amplitude=0.05, seed=seed)
            for _ in range(n_steps):
                sim.step()
            op[i, j] = float(sim.phi.std())
            Ld[i, j] = sim.domain_length()
            row.append(f"{op[i,j]:.2f}/{Ld[i,j]:.0f}")
        if verbose:
            print(f"  lambda={lam:.2f}: " + " ".join(row), flush=True)
    # classify
    op_thr = 0.18
    L_thr = L / 5.0
    reg = np.zeros_like(op)
    ordered = op > op_thr
    reg[ordered & (Ld < L_thr)] = 1                 # microphase
    reg[ordered & (Ld >= L_thr)] = 2                 # macrophase
    return dict(ks=ks, lams=lams, op=op, Ld=Ld, reg=reg, L=L)


# ===========================================================================
# Plane C:  (chi, N)  -- Flory-Huggins one-phase / two-phase (analytic)
# ===========================================================================

def plane_chi_N(n=320, Nmax=60, chi_max=3.5):
    """Two-phase gap (phi_R - phi_L from the spinodal) over (chi, N1) with the
    critical line chi_c(N)."""
    Ns = np.linspace(1.0, Nmax, n)
    chis = np.linspace(0.0, chi_max, n)
    CHI, NN = np.meshgrid(chis, Ns)
    # spinodal half-width: f'' = 1/(N phi) + 1/(1-phi) - 2 chi = 0 has two roots
    # when chi > chi_c(N); use the gap between them as the "two-phase strength".
    chic = 0.5 * (1.0 + 1.0 / np.sqrt(NN))**2
    gap = np.zeros_like(CHI)
    for idx in np.ndindex(CHI.shape):
        N1, chi = NN[idx], CHI[idx]
        if chi <= 0.5 * (1 + 1 / np.sqrt(N1))**2:
            continue
        # spinodal roots of 1/(N phi)+1/(1-phi)=2 chi
        coeffs = [2 * chi, -(2 * chi + 1 / N1 - 1), 1 / N1]   # in phi: A phi^2+..
        # 2chi phi^2 - (2chi + 1/N -1) phi + 1/N = 0  (from clearing denominators)
        r = np.roots([2 * chi, -(2 * chi - 1 + 1 / N1), 1 / N1])
        r = r[np.isreal(r)].real
        r = r[(r > 0) & (r < 1)]
        if len(r) == 2:
            gap[idx] = abs(r[1] - r[0])
    return dict(chis=chis, Ns=Ns, gap=gap, chic_line=(Ns, 0.5 * (1 + 1 / np.sqrt(Ns))**2))


# ===========================================================================
# Plane D:  (K, p)  -- rate-enhancement landscape (analytic)
# ===========================================================================

def plane_K_p(n=320, Kmax=300.0):
    """k_eff/k_out over (K, p) with K_A=K_B=K, plus the optimal ridge p*(K)."""
    K = np.logspace(0, np.log10(Kmax), n)
    p = np.linspace(1e-3, 1 - 1e-3, n)
    KK, PP = np.meshgrid(K, p)
    R = effective_rate_multiplier(KK, KK, PP)
    p_star = 1.0 / np.sqrt(K * K)                    # ~ 1/K interior optimum
    R_star = effective_rate_multiplier(K, K, np.clip(p_star, 1e-3, 1 - 1e-3))
    return dict(K=K, p=p, R=R, p_star=p_star, R_star=R_star)


# ===========================================================================
# Grand summary
# ===========================================================================

def plot_grand_summary(A, B, C, D, fname="param_space_grand.png"):
    set_style()
    fig, axes = plt.subplots(2, 2, figsize=(14.5, 12.0))
    (axA, axB), (axC, axD) = axes

    # ---- A: (phi_0, kappa) ----
    cmapR = ListedColormap(REG_COLORS[:3])
    axA.pcolormesh(A["phi0"], A["kappa"], A["reg"], cmap=cmapR, vmin=-0.5,
                   vmax=2.5, shading="auto")
    for x in (A["phi_s"], -A["phi_s"]):
        axA.axvline(x, color=INK, lw=1.6, ls="--")
    for x in (A["phi_b"], -A["phi_b"]):
        axA.axvline(x, color=INK, lw=1.4, ls=":")
    # finite-size cutoff lambda*=box
    axA.contour(A["phi0"], A["kappa"], A["lam_star"], levels=[A["box"]],
                colors=[ACCENT2], linewidths=2.0)
    axA.text(0.0, 7.0, "homogeneous\n(pattern $>$ box)", ha="center", fontsize=9,
             color=ACCENT2)
    axA.text(0.0, 2.0, "spinodal", ha="center", fontsize=12, color="white",
             fontweight="bold")
    axA.text(0.64, 4.0, "nucleation", ha="center", fontsize=10, color=INK,
             rotation=90)
    axA.text(0.90, 4.0, "homog.", ha="center", fontsize=9, color=INK, rotation=90)
    axA.set_xlabel(r"mean composition  $\phi_0$")
    axA.set_ylabel(r"gradient cost  $\kappa$  (box $L=%g$)" % A["box"])
    axA.set_title(r"(A) $(\phi_0,\kappa)$: spinodal / nucleation / homogeneous")
    # codim point where the spinodal edge (phi=0) meets the finite-size cutoff
    kap_cut = A["box"]**2 * 1.0 / (2 * np.pi)**2 / 2  # lambda*(0)=box -> kappa
    axA.plot(0.0, kap_cut, '*', color=GOLD, ms=20, mec=INK, mew=0.8, zorder=8)
    axA.annotate("codim-2:\nspinodal $\\to$ finite-size\ncutoff",
                 (0.0, kap_cut), (0.30, kap_cut - 1.8), fontsize=8,
                 color="#8a6a1a", arrowprops=dict(arrowstyle="->", color=GOLD))

    # ---- B: (k_fuel, lambda) ----
    cmapB = ListedColormap(REG_COLORS[:3])
    pcm = axB.pcolormesh(B["ks"], B["lams"], B["reg"], cmap=cmapB, vmin=-0.5,
                         vmax=2.5, shading="nearest")
    # overlay domain-length contours
    axB.contour(B["ks"], B["lams"], B["Ld"], levels=6, colors="k",
                linewidths=0.5, alpha=0.4)
    axB.set_xlabel(r"reaction rate  $k_{\rm fuel}$")
    axB.set_ylabel(r"activity  $\lambda$")
    axB.set_title(r"(B) $(k_{\rm fuel},\lambda)$: dissolution / microphase / "
                  "macrophase")
    # legend proxies
    from matplotlib.patches import Patch
    axB.legend(handles=[Patch(color=REG_COLORS[2], label="macrophase"),
                        Patch(color=REG_COLORS[1], label="microphase"),
                        Patch(color=REG_COLORS[0], label="dissolution")],
               loc="upper right", fontsize=8.5, framealpha=0.85)
    # codim-2 corner: where dissolution boundary meets micro/macro boundary
    axB.text(0.02, 0.12 * B["lams"].max(), "macrophase\n($\\lambda{=}0,k{=}0$)",
             fontsize=8, color="white")

    # ---- C: (chi, N) ----
    g = C["gap"].copy()
    g_masked = np.ma.masked_where(g <= 0, g)
    pcm2 = axC.pcolormesh(C["chis"], C["Ns"], g_masked, cmap="magma",
                          shading="auto")
    axC.plot(C["chic_line"][1], C["chic_line"][0], '-', color="cyan", lw=2.4,
             label=r"critical line $\chi_c(N)$")
    axC.set_xlabel(r"interaction  $\chi$")
    axC.set_ylabel(r"chain length  $N$")
    axC.set_title(r"(C) $(\chi,N)$: one-phase / two-phase (Flory--Huggins)")
    axC.text(2.6, 45, "two-phase", color="white", fontsize=11, ha="center")
    axC.text(0.55, 45, "one-phase", color=INK, fontsize=11, ha="center")
    axC.legend(loc="lower left", fontsize=9)
    cb2 = fig.colorbar(pcm2, ax=axC, fraction=0.046, pad=0.04)
    cb2.set_label(r"two-phase gap  $\phi_R-\phi_L$", fontsize=9)

    # ---- D: (K, p) ----
    pcm3 = axD.pcolormesh(D["K"], D["p"], D["R"], cmap="viridis",
                          norm=LogNorm(vmin=1.0, vmax=max(D["R"].max(), 2)),
                          shading="auto")
    axD.plot(D["K"], np.clip(D["p_star"], 1e-3, 1), '--', color="white",
             lw=2.0, label=r"optimal ridge  $p^*\!\approx\!1/K$")
    axD.set_xscale("log")
    axD.set_xlabel(r"partition coefficient  $K$")
    axD.set_ylabel(r"dense fraction  $p$")
    axD.set_title(r"(D) $(K,p)$: mass-conserving rate enhancement")
    axD.legend(loc="upper right", fontsize=9)
    cb3 = fig.colorbar(pcm3, ax=axD, fraction=0.046, pad=0.04)
    cb3.set_label(r"$k_{\rm eff}/k_{\rm out}$", fontsize=9)
    imax = np.unravel_index(np.argmax(D["R"]), D["R"].shape)
    axD.plot(D["K"][imax[1]], D["p"][imax[0]], '*', color=GOLD, ms=18,
             mec=INK, mew=0.7)

    fig.suptitle("The condensate model's parameter space: every accessible "
                 "regime in four control planes", fontsize=16, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    out = FIGURES / fname
    fig.savefig(out)
    plt.close(fig)
    print(f"  saved {out}")


if __name__ == "__main__":
    import sys
    quick = "--quick" in sys.argv

    print("Parameter-space exploration")
    print("=" * 40)
    print("\n[A] (phi_0, kappa) regime map (analytic)")
    A = plane_phi0_kappa()
    print("[C] (chi, N) Flory-Huggins (analytic)")
    C = plane_chi_N(n=200 if quick else 320)
    print("[D] (K, p) rate enhancement (analytic)")
    D = plane_K_p(n=200 if quick else 320)
    print("[B] (k_fuel, lambda) active-regime sweep (simulation)")
    B = plane_fuel_activity(nk=7, nl=5, N=48, T=160.0) if quick \
        else plane_fuel_activity(nk=11, nl=9, N=64, T=260.0)

    print("\n[*] Grand summary figure")
    plot_grand_summary(A, B, C, D)
    print("\nDone.")
