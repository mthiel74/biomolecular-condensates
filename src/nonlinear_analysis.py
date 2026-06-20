#!/usr/bin/env python3
r"""Nonlinear-dynamics analysis of the Cahn–Hilliard + reaction model.

This module treats the condensate model already implemented in this repo
(`flory_huggins.py`, `cahn_hilliard.py`, `condensate_kinetics.py`) with the
standard toolkit of nonlinear dynamics:

1.  Linear stability of the conserved Cahn–Hilliard dynamics
        ∂φ/∂t = M ∇²( f'(φ) − κ ∇²φ )
    A perturbation φ = φ₀ + δ·exp(i k·x + ω t) grows at rate
        ω(k) = −M k² ( f''(φ₀) + κ k² ).
    The band 0 < k < k_c = sqrt(−f''(φ₀)/κ) is unstable whenever
    f''(φ₀) < 0 (the spinodal).  The fastest-growing mode and its rate are
        k*   = sqrt( −f''(φ₀) / (2κ) ),     λ* = 2π/k*,
        ω*   = M f''(φ₀)² / (4κ).
    We instantiate f'' both for the clean symmetric double well
    f(φ)=−a/2 φ²+b/4 φ⁴  (the canonical pedagogical limit) and for the
    Flory–Huggins free energy that the repo's solver actually integrates.

2.  Bifurcation / phase diagram via the common-tangent (Maxwell)
    construction:  the dilute/dense coexistence fractions as the control
    parameter is swept.  For the double well this is the textbook
    pitchfork-shaped binodal φ_b = ±sqrt(a/b) enclosing the spinodal
    φ_s = ±sqrt(a/3b);  for Flory–Huggins it is the χ–φ phase diagram.

3.  Coarsening:  the structure-factor length scale L(t) extracted from the
    repo's own solver, tested against the Lifshitz–Slyozov–Wagner t^{1/3}.

4.  Energy landscape:  the Lyapunov functional F[φ] and its bulk/gradient
    split decreasing monotonically along the trajectory.

5.  Phase portrait:  the trajectory in (dense-phase area fraction, F) space
    relaxing onto the equilibrium fixed point.

6.  Reaction–diffusion coupling:  a mass-conserving "rapid-partition" model
    for the effective bimolecular rate, k_eff/k_out, as a function of the
    partition coefficients (K_A, K_B) and the dense area fraction p.

All figures are written to ``figures/`` at PRL quality.

Run as a script to regenerate every figure:  ``python3 src/nonlinear_analysis.py``
"""

from __future__ import annotations

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pathlib import Path

# Repo modules ---------------------------------------------------------------
from flory_huggins import d2f_dphi2, critical_point
from cahn_hilliard import CahnHilliard2D

# ---------------------------------------------------------------------------
# A correctly-coarsening Model-B solver (symmetric double well).
#
# NOTE (honest caveat): the repo's Flory–Huggins solver uses an Eyre
# stabilisation constant A=50 — forced by the log singularity of f(φ) as
# φ→0.  Such heavy stabilisation injects artificial dissipation that *pins
# interfaces and arrests coarsening*: in a large box L(t) freezes after the
# initial instability instead of following the LSW t^{1/3} law.  The
# canonical demonstration of Lifshitz–Slyozov–Wagner coarsening therefore
# uses the symmetric double well f(φ)=−a/2 φ²+b/4 φ⁴, which has no
# singularity, needs only modest stabilisation (A≈2.5), and reproduces
# t^{1/3} cleanly.  We use it for the dynamical studies (coarsening, energy
# landscape, phase portrait) and keep Flory–Huggins for the thermodynamics
# (dispersion, phase diagram).
# ---------------------------------------------------------------------------


class DoubleWellCH2D:
    r"""Semi-implicit spectral Model-B solver for the double-well free energy.

    f(φ) = −a/2 φ² + b/4 φ⁴,   μ = f'(φ) = −a φ + b φ³,
    ∂φ/∂t = M ∇²(μ − κ ∇²φ).

    Eyre splitting μ = (μ − Aφ) + Aφ with A ≳ max f'' treated implicitly:
        φ̂^{n+1} = (φ̂ − dt M k² ĝ) / (1 + dt M k² (A + κ k²)),  g = μ − Aφ.
    """

    def __init__(self, L=128.0, N=256, a=1.0, b=1.0, kappa=1.0,
                 mobility=1.0, dt=0.05, A=2.5):
        self.L, self.N = L, N
        self.a, self.b, self.kappa = a, b, kappa
        self.M, self.dt, self.A = mobility, dt, A
        self.dx = L / N
        k = 2 * np.pi * np.fft.fftfreq(N, d=self.dx)
        KX, KY = np.meshgrid(k, k)
        self.k2 = KX**2 + KY**2
        self._denom = 1.0 + dt * mobility * self.k2 * (A + kappa * self.k2)
        self.phi = None
        self.time = 0.0

    def initialize(self, phi_mean=0.0, noise_amplitude=0.02, seed=42):
        rng = np.random.default_rng(seed)
        self.phi = phi_mean + noise_amplitude * rng.standard_normal((self.N, self.N))
        self.time = 0.0

    def step(self):
        mu = -self.a * self.phi + self.b * self.phi**3
        g = mu - self.A * self.phi
        phi_hat = np.fft.fft2(self.phi)
        g_hat = np.fft.fft2(g)
        self.phi = np.real(np.fft.ifft2(
            (phi_hat - self.dt * self.M * self.k2 * g_hat) / self._denom))
        self.time += self.dt

FIGURES = Path(__file__).resolve().parent.parent / "figures"
FIGURES.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Shared aesthetic — a single consistent "PRL" style for every figure.
# ---------------------------------------------------------------------------

# A perceptually-uniform accent ramp reused across panels.
SEQ_CMAP = "magma"
DIV_CMAP = "RdBu_r"
ACCENT = "#c81e5b"      # dense / dramatic
ACCENT2 = "#1f6fb2"     # dilute / cool
GOLD = "#e8a33d"
INK = "#1a1a1a"


def set_style():
    mpl.rcParams.update({
        "figure.dpi": 120,
        "savefig.dpi": 240,
        "savefig.bbox": "tight",
        "font.family": "serif",
        "font.serif": ["DejaVu Serif", "Times New Roman", "Computer Modern Roman"],
        "mathtext.fontset": "cm",
        "axes.titlesize": 14,
        "axes.labelsize": 13,
        "axes.linewidth": 1.0,
        "axes.edgecolor": "#333333",
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "legend.frameon": False,
        "legend.fontsize": 10,
        "lines.linewidth": 2.0,
        "axes.grid": True,
        "grid.alpha": 0.18,
        "grid.linewidth": 0.6,
    })


# ===========================================================================
# Robust Flory–Huggins binodal
#
# The repo's binodal_curve() uses a 2-D fsolve that fails for the strongly
# asymmetric case (N1=100, N2=1): it returns only ~6 points before the dilute
# branch φ_L underflows.  Here we use the structure of the problem instead:
#   • dense branch  φ_R(χ):  root of the osmotic pressure Π(φ,χ)=0  — exact in
#     the limit φ_L→0 (dense phase coexisting with near-pure solvent), which
#     holds for all χ more than a hair above χ_c here;
#   • dilute branch φ_L(χ):  μ(φ_L,χ)=μ(φ_R,χ) solved in log φ.
# Verified: μ(φ_L)=μ(φ_R) to 4 decimals over the whole range.
# ===========================================================================

def robust_binodal(N1=100, N2=1, chi_max_factor=3.0, n=160):
    """Coexistence curve over χ ∈ (χ_c, chi_max_factor·χ_c)."""
    from scipy.optimize import brentq
    from flory_huggins import osmotic_pressure, chemical_potential
    phi_c, chi_c = critical_point(N1, N2)
    chis = np.linspace(chi_c * 1.01, chi_c * chi_max_factor, n)
    out_chi, out_L, out_R = [], [], []
    for chi in chis:
        try:
            phiR = brentq(lambda p: osmotic_pressure(p, chi, N1, N2),
                          phi_c + 1e-6, 1 - 1e-10)
            mu_c = chemical_potential(phiR, chi, N1, N2)
            f = lambda lp: chemical_potential(np.exp(lp), chi, N1, N2) - mu_c
            lo, hi = np.log(1e-300), np.log(phi_c - 1e-6)
            phiL = np.exp(brentq(f, lo, hi)) if f(lo) * f(hi) < 0 else np.nan
        except (ValueError, RuntimeError):
            continue
        out_chi.append(chi); out_L.append(phiL); out_R.append(phiR)
    return np.array(out_chi), np.array(out_L), np.array(out_R)


# ===========================================================================
# 1. Linear stability — dispersion relation ω(k)
# ===========================================================================

def double_well_fpp(phi, a=1.0, b=1.0):
    r"""Second derivative of the symmetric double well f=−a/2 φ²+b/4 φ⁴."""
    return -a + 3 * b * phi**2


def dispersion(k, fpp, M=1.0, kappa=0.5):
    r"""Cahn–Hilliard growth rate ω(k) = −M k² (f'' + κ k²)."""
    return -M * k**2 * (fpp + kappa * k**2)


def fastest_mode(fpp, kappa=0.5, M=1.0):
    r"""Return (k*, ω*) of the fastest-growing mode, or (nan, nan) if stable."""
    if fpp >= 0:
        return np.nan, np.nan
    kstar = np.sqrt(-fpp / (2 * kappa))
    omega = M * fpp**2 / (4 * kappa)
    return kstar, omega


def plot_dispersion_relation():
    """ω(k) for the double well and for Flory–Huggins, with k* marked."""
    set_style()
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.5, 5.4))
    k = np.linspace(0, 3.0, 600)
    kappa, M = 0.5, 1.0

    # --- Panel A: symmetric double well, vary mean composition φ₀ ----------
    phi0_vals = np.linspace(0.0, 0.62, 6)
    cmap = mpl.colormaps[SEQ_CMAP]
    for i, phi0 in enumerate(phi0_vals):
        fpp = double_well_fpp(phi0, a=1.0, b=1.0)
        w = dispersion(k, fpp, M, kappa)
        c = cmap(0.12 + 0.72 * i / (len(phi0_vals) - 1))
        axA.plot(k, w, color=c, lw=2.2,
                 label=fr"$\phi_0={phi0:.2f}$  ($f''={fpp:+.2f}$)")
        kstar, wstar = fastest_mode(fpp, kappa, M)
        if np.isfinite(kstar):
            axA.plot(kstar, wstar, 'o', color=c, ms=6, zorder=5)
    axA.axhline(0, color=INK, lw=0.8)
    # spinodal edge for the double well: |phi0| = sqrt(a/3b)
    axA.axvline(0, color="none")
    axA.set_title(r"(a) Symmetric double well  $f=-\frac{a}{2}\phi^2+\frac{b}{4}\phi^4$")
    axA.set_xlabel(r"Wavenumber  $k$")
    axA.set_ylabel(r"Growth rate  $\omega(k)$")
    axA.set_xlim(0, 3)
    axA.set_ylim(-1.2, 0.4)
    axA.legend(loc="lower left", fontsize=8.5)
    axA.text(0.02, 0.95, "unstable band\n" r"$0<k<k_c$", transform=axA.transAxes,
             va="top", fontsize=9, color=ACCENT)

    # --- Panel B: Flory–Huggins (the model the solver integrates) ----------
    chi, N1, N2 = 1.5, 100, 1
    phi0_vals = [0.06, 0.10, 0.20, 0.30, 0.45, 0.60]
    for i, phi0 in enumerate(phi0_vals):
        fpp = d2f_dphi2(phi0, chi, N1, N2)
        w = dispersion(k, fpp, M, kappa)
        c = cmap(0.12 + 0.72 * i / (len(phi0_vals) - 1))
        axB.plot(k, w, color=c, lw=2.2,
                 label=fr"$\phi_0={phi0:.2f}$  ($f''={fpp:+.2f}$)")
        kstar, wstar = fastest_mode(fpp, kappa, M)
        if np.isfinite(kstar):
            axB.plot(kstar, wstar, 'o', color=c, ms=6, zorder=5)
    axB.axhline(0, color=INK, lw=0.8)
    axB.set_title(fr"(b) Flory–Huggins  ($\chi={chi}$, $N_1={N1}$, $N_2={N2}$)")
    axB.set_xlabel(r"Wavenumber  $k$")
    axB.set_ylabel(r"Growth rate  $\omega(k)$")
    axB.set_xlim(0, 3)
    axB.set_ylim(-1.6, 1.6)
    axB.legend(loc="lower left", fontsize=8.5)

    fig.suptitle("Linear stability of Cahn–Hilliard dynamics: dispersion relation "
                 r"$\omega(k)=-Mk^2\,(f''(\phi_0)+\kappa k^2)$",
                 fontsize=14.5, y=1.02)
    fig.tight_layout()
    out = FIGURES / "nld_dispersion_relation.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"  saved {out}")


def plot_fastest_mode_map():
    r"""k*(φ₀, κ) heatmap with the spinodal boundary, for Flory–Huggins."""
    set_style()
    chi, N1, N2 = 1.5, 100, 1
    phi = np.linspace(0.005, 0.75, 400)
    kappa = np.linspace(0.05, 2.0, 400)
    PHI, KAP = np.meshgrid(phi, kappa)
    FPP = d2f_dphi2(PHI, chi, N1, N2)

    # fastest-growing wavelength λ* = 2π/k*, defined only where f''<0
    with np.errstate(invalid="ignore", divide="ignore"):
        KSTAR = np.sqrt(-FPP / (2 * KAP))
        LAMBDA = 2 * np.pi / KSTAR
    LAMBDA = np.ma.masked_where(FPP >= 0, LAMBDA)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.4))

    # Panel A: λ*(φ₀, κ)
    pcm = ax1.pcolormesh(PHI, KAP, LAMBDA, cmap="viridis",
                         norm=LogNorm(vmin=2, vmax=40), shading="auto")
    cb = fig.colorbar(pcm, ax=ax1)
    cb.set_label(r"Fastest-growing wavelength  $\lambda^*=2\pi/k^*$")
    # spinodal boundary f''=0
    ax1.contour(PHI, KAP, FPP, levels=[0], colors="white", linewidths=2.0)
    ax1.set_title(r"(a) Selected length scale  $\lambda^*(\phi_0,\kappa)$")
    ax1.set_xlabel(r"Mean composition  $\phi_0$")
    ax1.set_ylabel(r"Gradient coefficient  $\kappa$")
    ax1.text(0.62, 1.5, "stable\n" r"($f''>0$)", color="white", fontsize=11,
             ha="center", va="center")

    # Panel B: spinodal region in (φ₀, κ) — actually f''(φ₀) is κ-independent,
    # so the *unstable band exists* for φ between the spinodal roots regardless
    # of κ; κ sets the band width k_c = sqrt(-f''/κ).  Plot k_c.
    with np.errstate(invalid="ignore"):
        KC = np.sqrt(-FPP / KAP)
    KC = np.ma.masked_where(FPP >= 0, KC)
    pcm2 = ax2.pcolormesh(PHI, KAP, KC, cmap="magma", shading="auto")
    cb2 = fig.colorbar(pcm2, ax=ax2)
    cb2.set_label(r"Cutoff wavenumber  $k_c=\sqrt{-f''/\kappa}$")
    ax2.contour(PHI, KAP, FPP, levels=[0], colors="cyan", linewidths=2.0)
    # mark the two spinodal compositions
    from flory_huggins import _spinodal_roots
    roots = _spinodal_roots(chi, N1, N2)
    if roots:
        for r in roots:
            ax2.axvline(r, color="cyan", ls=":", lw=1.2)
        ax2.text(0.5 * (roots[0] + roots[1]), 0.2, "spinodal band",
                 color="cyan", ha="center", fontsize=10)
    ax2.set_title(r"(b) Unstable band width  $k_c(\phi_0,\kappa)$")
    ax2.set_xlabel(r"Mean composition  $\phi_0$")
    ax2.set_ylabel(r"Gradient coefficient  $\kappa$")

    fig.suptitle(fr"Mode selection across $(\phi_0,\kappa)$  —  Flory–Huggins, "
                 fr"$\chi={chi}$, $N_1={N1}$", fontsize=14.5, y=1.02)
    fig.tight_layout()
    out = FIGURES / "nld_mode_selection_map.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"  saved {out}")


# ===========================================================================
# 2. Bifurcation diagram / phase diagram
# ===========================================================================

def plot_bifurcation_diagram():
    """Coexistence (binodal) and stability (spinodal) branches.

    Panel (a): canonical symmetric double well, control parameter a (well
    depth).  Binodal φ=±sqrt(a/b), spinodal φ=±sqrt(a/3b) — a supercritical
    pitchfork in the coexistence amplitude.
    Panel (b): the repo's Flory–Huggins free energy, control parameter χ,
    binodal from the common-tangent solver in flory_huggins.py.
    """
    set_style()
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.5, 5.6))

    # --- Panel A: double well ---------------------------------------------
    b = 1.0
    a = np.linspace(0, 1.5, 400)
    phi_bin = np.sqrt(np.clip(a / b, 0, None))
    phi_spin = np.sqrt(np.clip(a / (3 * b), 0, None))

    axA.plot(a, phi_bin, color=ACCENT2, lw=2.6, label="binodal (coexistence)")
    axA.plot(a, -phi_bin, color=ACCENT2, lw=2.6)
    axA.plot(a, phi_spin, color=ACCENT, lw=2.2, ls="--", label="spinodal")
    axA.plot(a, -phi_spin, color=ACCENT, lw=2.2, ls="--")
    axA.fill_between(a, -phi_spin, phi_spin, color=ACCENT, alpha=0.10)
    axA.fill_between(a, phi_spin, phi_bin, color=GOLD, alpha=0.18)
    axA.fill_between(a, -phi_bin, -phi_spin, color=GOLD, alpha=0.18)
    axA.axhline(0, color=INK, lw=0.8)
    axA.plot(0, 0, 'o', color=INK, ms=8, zorder=6)
    axA.annotate("critical point\n(2nd-order)", xy=(0, 0), xytext=(0.45, 0.18),
                 fontsize=9, arrowprops=dict(arrowstyle="->", color=INK))
    axA.text(1.25, 0.05, "unstable\n(spinodal)", color=ACCENT, ha="center",
             fontsize=9)
    axA.text(1.3, 0.78, "metastable", color="#9a7414", ha="center", fontsize=9)
    axA.text(1.3, 1.18, "stable single phase",
             color=ACCENT2, ha="center", fontsize=9)
    axA.set_title(r"(a) Double well: pitchfork in coexistence amplitude")
    axA.set_xlabel(r"Control parameter  $a$  (well depth $\propto\,\chi-\chi_c$)")
    axA.set_ylabel(r"Order parameter  $\phi-\phi_c$")
    axA.set_xlim(0, 1.5)
    axA.set_ylim(-1.35, 1.35)
    axA.legend(loc="upper left")

    # --- Panel B: Flory–Huggins phase diagram (χ vs φ) --------------------
    N1, N2 = 100, 1
    phi_c, chi_c = critical_point(N1, N2)
    from flory_huggins import spinodal_curve
    phi_sp, chi_sp = spinodal_curve(N1, N2, n_points=400)
    chi_bn, phL, phR = robust_binodal(N1, N2, chi_max_factor=3.0, n=200)

    if len(chi_bn):
        axB.plot(phL, chi_bn, color=ACCENT2, lw=2.6, label="binodal")
        axB.plot(phR, chi_bn, color=ACCENT2, lw=2.6)
        axB.fill_betweenx(chi_bn, phL, phR, color=ACCENT2, alpha=0.07)
    mid = len(phi_sp) // 2
    axB.plot(phi_sp[:mid], chi_sp[:mid], color=ACCENT, ls="--", lw=2.2,
             label="spinodal")
    axB.plot(phi_sp[mid:], chi_sp[mid:], color=ACCENT, ls="--", lw=2.2)
    axB.fill_betweenx(chi_sp[:mid], phi_sp[:mid], phi_sp[mid:][::-1],
                      color=ACCENT, alpha=0.08)
    axB.plot(phi_c, chi_c, 'o', color=INK, ms=8, zorder=6)
    axB.annotate(fr"critical point" "\n" fr"$(\phi_c={phi_c:.3f},\ \chi_c={chi_c:.3f})$",
                 xy=(phi_c, chi_c), xytext=(phi_c + 0.18, chi_c + 0.25),
                 fontsize=9, arrowprops=dict(arrowstyle="->", color=INK))
    axB.text(0.35, 1.4, "two-phase\n(dilute + dense)", color=ACCENT2,
             ha="center", fontsize=9)
    axB.text(0.35, 0.45, "one phase (mixed)", color="green", ha="center",
             fontsize=9)
    # the simulation operating point
    axB.plot(0.30, 1.5, '*', color=GOLD, ms=16, mec=INK, mew=0.6, zorder=7,
             label=r"solver: $\phi_0=0.3,\ \chi=1.5$")
    axB.set_title(fr"(b) Flory–Huggins phase diagram ($N_1={N1}$)")
    axB.set_xlabel(r"Volume fraction  $\phi$")
    axB.set_ylabel(r"Interaction  $\chi\ (\propto 1/T)$")
    axB.set_xlim(0, 0.8)
    axB.set_ylim(0, chi_c * 3)
    axB.legend(loc="upper right")

    fig.suptitle("Bifurcation structure: common-tangent coexistence vs. "
                 "spinodal instability", fontsize=14.5, y=1.02)
    fig.tight_layout()
    out = FIGURES / "nld_bifurcation_diagram.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"  saved {out}")


# ===========================================================================
# Shared CH diagnostics — run the solver ONCE, record every observable.
# ===========================================================================

def _radial_structure_factor(phi, k2, dx):
    """Return (k_centers, S_radial) for the fluctuation field."""
    f = phi - phi.mean()
    S = np.abs(np.fft.fft2(f))**2
    kmag = np.sqrt(k2)
    kmax = kmag.max() / 2
    bins = np.linspace(0, kmax, 60)
    Sr = np.zeros(len(bins) - 1)
    for i in range(len(bins) - 1):
        m = (kmag >= bins[i]) & (kmag < bins[i + 1])
        if m.any():
            Sr[i] = S[m].mean()
    kc = 0.5 * (bins[:-1] + bins[1:])
    return kc, Sr


def _domain_length_first_moment(phi, k2):
    r"""Coarsening length L = 2π / <k>, first moment of S(k)."""
    f = phi - phi.mean()
    S = np.abs(np.fft.fft2(f))**2
    kmag = np.sqrt(k2)
    sel = kmag > 0
    k1 = (kmag[sel] * S[sel]).sum() / S[sel].sum()
    return 2 * np.pi / k1 if k1 > 0 else np.nan


def _energy_components(sim):
    """Bulk and gradient contributions to F[φ] for the double-well solver."""
    phi = sim.phi
    f_bulk = -0.5 * sim.a * phi**2 + 0.25 * sim.b * phi**4
    gx = (np.roll(phi, -1, 1) - np.roll(phi, 1, 1)) / (2 * sim.dx)
    gy = (np.roll(phi, -1, 0) - np.roll(phi, 1, 0)) / (2 * sim.dx)
    f_grad = 0.5 * sim.kappa * (gx**2 + gy**2)
    F_bulk = f_bulk.sum() * sim.dx**2
    F_grad = f_grad.sum() * sim.dx**2
    return F_bulk, F_grad


def run_diagnostics(total_time=1500.0, N=256, L=128.0, phi_mean=0.0,
                    kappa=1.0, dt=0.05, n_record=150):
    """Single double-well Model-B run: records L(t), F(t) split, dense fraction.

    The symmetric double well coarsens correctly (see DoubleWellCH2D note),
    so this single trajectory feeds the coarsening, energy-landscape and
    phase-portrait figures.
    """
    print(f"  Running double-well Model-B solver  (N={N}, L={L}, T={total_time}, "
          f"~{int(total_time/dt)} steps)...")
    sim = DoubleWellCH2D(L=L, N=N, kappa=kappa, dt=dt)
    sim.initialize(phi_mean=phi_mean)

    # "dense phase" = region settled into the +1 well (φ > 0.5)
    thresh = 0.5

    n_total = int(total_time / dt)
    rec_every = max(1, n_total // n_record)

    rec = {k: [] for k in ("t", "L", "F", "F_bulk", "F_grad", "dense_frac")}
    snaps = []
    snap_times = [0.0, 2.0, 20.0, 150.0, 600.0, total_time]

    def record():
        Fb, Fg = _energy_components(sim)
        rec["t"].append(sim.time)
        rec["L"].append(_domain_length_first_moment(sim.phi, sim.k2))
        rec["F_bulk"].append(Fb)
        rec["F_grad"].append(Fg)
        rec["F"].append(Fb + Fg)
        rec["dense_frac"].append(float((sim.phi > thresh).mean()))

    record()
    snaps.append((0.0, sim.phi.copy()))
    next_snap = 1
    for i in range(1, n_total + 1):
        sim.step()
        if i % rec_every == 0:
            record()
        if next_snap < len(snap_times) and sim.time >= snap_times[next_snap]:
            snaps.append((sim.time, sim.phi.copy()))
            next_snap += 1
    # always capture the final state (float drift can miss the last threshold)
    if snaps[-1][0] < sim.time - 0.5 * dt:
        snaps.append((sim.time, sim.phi.copy()))

    for k in rec:
        rec[k] = np.array(rec[k])
    rec["thresh"] = thresh
    return sim, rec, snaps


def fh_frozen_coarsening(total_time=150.0, N=128, L=40.0, chi=1.5,
                         phi_mean=0.3, kappa=0.5, dt=0.005):
    """Run the repo's Flory–Huggins solver and return (t, L) — demonstrates
    the coarsening arrest caused by heavy Eyre stabilisation (A=50)."""
    sim = CahnHilliard2D(L=L, N=N, chi=chi, kappa=kappa, dt=dt)
    sim.initialize(phi_mean=phi_mean)
    n_total = int(total_time / dt)
    rec_every = max(1, n_total // 80)
    t, Lt = [], []
    for i in range(1, n_total + 1):
        sim.step()
        if i % rec_every == 0:
            t.append(sim.time)
            Lt.append(_domain_length_first_moment(sim.phi, sim.k2))
    return np.array(t), np.array(Lt)


# ===========================================================================
# 3. Coarsening dynamics — L(t) ~ t^{1/3}
# ===========================================================================

def plot_coarsening(rec, fh_curve=None):
    """Two-panel LSW coarsening figure.

    (a) log-log L(t) with the t^{1/3} guide; the frozen Flory–Huggins solver
        is overlaid to show coarsening arrest.
    (b) the rigorous LSW statement L³ = L₀³ + K t as a straight line.
    """
    set_style()
    t, L = rec["t"], rec["L"]
    m = t > 5.0  # skip the early linear-instability transient
    t, L = t[m], L[m]

    # apparent log-log slope on the late window (offset-contaminated)
    late = t > 0.2 * t.max()
    alpha = np.polyfit(np.log(t[late]), np.log(L[late]), 1)[0]
    # rigorous LSW fit: L³ linear in t
    K, L03 = np.polyfit(t, L**3, 1)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.0, 5.4))

    # --- Panel A: log-log -------------------------------------------------
    axA.loglog(t, L, 'o', color=ACCENT, ms=5, mec="white", mew=0.5,
               label="double-well Model B")
    tt = np.linspace(t.min(), t.max(), 100)
    cref = L[len(L)//2] / t[len(t)//2]**(1/3)
    axA.loglog(tt, cref * tt**(1/3), '--', color=INK, lw=2.2,
               label=r"LSW  $L\sim t^{1/3}$")
    if fh_curve is not None:
        tf, Lf = fh_curve
        good = tf > 2
        axA.loglog(tf[good], Lf[good], 's-', color=ACCENT2, ms=3, lw=1.2,
                   alpha=0.8, label="Flory–Huggins solver (frozen)")
    axA.set_xlabel(r"Time  $t$")
    axA.set_ylabel(r"Domain length  $L(t)=2\pi/\langle k\rangle$")
    axA.set_title("(a) Coarsening: power-law growth")
    axA.legend(loc="upper left", fontsize=9)

    # --- Panel B: L³ vs t (the rigorous LSW law) --------------------------
    axB.plot(t, L**3, 'o', color=ACCENT, ms=5, mec="white", mew=0.5,
             label=r"$L^3(t)$")
    axB.plot(t, K * t + L03, '-', color=INK, lw=2.2,
             label=fr"$L^3=L_0^3+Kt,\ K={K:.1f}$")
    axB.set_xlabel(r"Time  $t$")
    axB.set_ylabel(r"$L^3(t)$")
    axB.set_title(r"(b) Rigorous LSW:  $L^3$ linear in $t$")
    axB.legend(loc="upper left")
    # R² of the linear fit
    pred = K * t + L03
    ss_res = np.sum((L**3 - pred)**2)
    ss_tot = np.sum((L**3 - np.mean(L**3))**2)
    r2 = 1 - ss_res / ss_tot
    axB.text(0.96, 0.08,
             fr"$R^2={r2:.4f}$" "\n" fr"log-log slope $={alpha:.2f}$" "\n"
             r"($\to 1/3$ asymptotically)",
             transform=axB.transAxes, ha="right", va="bottom", fontsize=10,
             bbox=dict(boxstyle="round", fc="white", ec="#cccccc", alpha=0.9))

    fig.suptitle("Lifshitz–Slyozov–Wagner coarsening of the conserved order "
                 "parameter", fontsize=14.5, y=1.02)
    fig.tight_layout()
    out = FIGURES / "nld_coarsening_scaling.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"  saved {out}  (L^3 slope K={K:.2f}, R^2={r2:.4f}, "
          f"log-log alpha={alpha:.3f})")
    return alpha


# ===========================================================================
# 4. Energy landscape — F(t) and its decomposition
# ===========================================================================

def plot_energy_landscape(rec):
    set_style()
    t = rec["t"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.0, 5.2))

    ax1.plot(t, rec["F"], color=INK, lw=2.6, label=r"total  $F$")
    ax1.plot(t, rec["F_bulk"], color=ACCENT2, lw=2.2, label=r"bulk  $\int f(\phi)$")
    ax1.plot(t, rec["F_grad"], color=ACCENT, lw=2.2,
             label=r"gradient  $\int \frac{\kappa}{2}|\nabla\phi|^2$")
    ax1.set_xlabel(r"Time  $t$")
    ax1.set_ylabel(r"Free energy  $F[\phi]$")
    ax1.set_title("(a) Lyapunov functional decreases monotonically")
    ax1.legend(loc="center right")

    # monotonicity check: dF/dt should be <= 0
    dF = np.gradient(rec["F"], t)
    ax2.plot(t, dF, color=ACCENT, lw=2.0)
    ax2.axhline(0, color=INK, lw=0.9, ls="--")
    ax2.fill_between(t, dF, 0, where=dF <= 0, color=ACCENT2, alpha=0.15)
    ax2.set_xlabel(r"Time  $t$")
    ax2.set_ylabel(r"$dF/dt$")
    ax2.set_title(r"(b) Dissipation rate  $dF/dt\leq 0$")
    frac_neg = float((dF <= 1e-9).mean())
    ax2.text(0.96, 0.08, fr"$dF/dt\leq0$ for {100*frac_neg:.0f}% of steps",
             transform=ax2.transAxes, ha="right", va="bottom", fontsize=10,
             bbox=dict(boxstyle="round", fc="white", ec="#cccccc", alpha=0.9))

    fig.suptitle("Energy landscape of Cahn–Hilliard relaxation",
                 fontsize=14.5, y=1.01)
    fig.tight_layout()
    out = FIGURES / "nld_energy_landscape.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"  saved {out}")


# ===========================================================================
# 5. Phase portrait — (dense fraction, F)
# ===========================================================================

def plot_phase_portrait(rec, snaps):
    set_style()
    fig, ax = plt.subplots(figsize=(8.2, 6.2))

    x = rec["dense_frac"]
    y = rec["F"]
    t = rec["t"]

    pts = np.array([x, y]).T.reshape(-1, 1, 2)
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
    from matplotlib.collections import LineCollection
    lc = LineCollection(segs, cmap="magma",
                        norm=mpl.colors.Normalize(t.min(), t.max()))
    lc.set_array(t[:-1])
    lc.set_linewidth(3.0)
    ax.add_collection(lc)
    cb = fig.colorbar(lc, ax=ax)
    cb.set_label(r"Time  $t$")

    ax.plot(x[0], y[0], 'o', color=ACCENT2, ms=11, mec="white", zorder=6,
            label="initial (mixed + noise)")
    ax.plot(x[-1], y[-1], '*', color=ACCENT, ms=20, mec=INK, mew=0.6, zorder=6,
            label="fixed point (phase-separated)")

    ax.set_xlabel(r"Dense-phase area fraction  $p(t)$")
    ax.set_ylabel(r"Free energy  $F[\phi]$")
    ax.set_title("Phase portrait: relaxation onto the equilibrium manifold")
    ax.legend(loc="upper right")
    ax.margins(0.08)
    fig.tight_layout()
    out = FIGURES / "nld_phase_portrait.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"  saved {out}")


# ===========================================================================
# 6. Reaction–diffusion coupling — mass-conserving rapid-partition rate
# ===========================================================================

def effective_rate_multiplier(K_A, K_B, p):
    r"""Mass-conserving bimolecular rate enhancement for A+B→C.

    With dense area fraction p and rapid partitioning c_in = K c_out in each
    phase, conservation of *total* A and B fixes the dilute concentrations.
    The volume-averaged reaction rate, relative to the same total amounts in a
    well-mixed system, is

        k_eff/k_out =  [p K_A K_B + (1-p)]
                       ----------------------------------------
                       (1 + p(K_A-1)) (1 + p(K_B-1))

    Limits:  p→0 and p→1 both give 1 (no net gain); the maximum is interior,
    ≈ sqrt(K_A K_B)/4 near p ≈ 1/sqrt(K_A K_B).  This is the honest
    correction to the naive "K_A·K_B" enhancement, which is the open-system
    (buffered-reservoir) limit where the dilute phase never depletes.
    """
    num = p * K_A * K_B + (1 - p)
    den = (1 + p * (K_A - 1)) * (1 + p * (K_B - 1))
    return num / den


def plot_reaction_coupling():
    set_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.6))

    # --- Panel A: k_eff/k_out over (K_A, K_B) at fixed dense fraction ------
    p = 0.05
    K = np.logspace(0, 3, 300)
    KA, KB = np.meshgrid(K, K)
    E = effective_rate_multiplier(KA, KB, p)
    pcm = ax1.pcolormesh(KA, KB, E, cmap=SEQ_CMAP,
                         norm=LogNorm(vmin=1, vmax=E.max()), shading="auto")
    ax1.set_xscale("log"); ax1.set_yscale("log")
    cb = fig.colorbar(pcm, ax=ax1)
    cb.set_label(r"$k_{\mathrm{eff}}/k_{\mathrm{out}}$")
    cs = ax1.contour(KA, KB, E, levels=[2, 5, 10, 20], colors="white",
                     linewidths=1.0)
    ax1.clabel(cs, fmt="%.0f", fontsize=8)
    ax1.set_xlabel(r"Partition coefficient  $K_A$")
    ax1.set_ylabel(r"Partition coefficient  $K_B$")
    ax1.set_title(fr"(a) Enhancement over $(K_A,K_B)$  at  $p={p}$")

    # --- Panel B: k_eff/k_out vs dense fraction p, several K --------------
    pp = np.linspace(1e-4, 1, 600)
    cmap = mpl.colormaps["viridis"]
    Kvals = [10, 30, 100, 300, 1000]
    for i, Kv in enumerate(Kvals):
        E = effective_rate_multiplier(Kv, Kv, pp)
        c = cmap(0.1 + 0.8 * i / (len(Kvals) - 1))
        ax2.plot(pp, E, color=c, lw=2.4, label=fr"$K_A=K_B={Kv}$")
        # mark optimum
        imax = np.argmax(E)
        ax2.plot(pp[imax], E[imax], 'o', color=c, ms=6, zorder=5)
    ax2.axhline(1, color=INK, lw=0.9, ls=":")
    ax2.set_yscale("log")
    ax2.set_xlabel(r"Dense-phase area fraction  $p$")
    ax2.set_ylabel(r"$k_{\mathrm{eff}}/k_{\mathrm{out}}$")
    ax2.set_title(r"(b) Optimal partitioning at intermediate $p^\ast\approx 1/K$")
    ax2.legend(loc="upper right", fontsize=9)
    ax2.text(0.5, 0.05,
             "naive $K^2$ unreachable:\nbulk depletion caps gain at $\\sim K/4$",
             transform=ax2.transAxes, ha="center", va="bottom", fontsize=9,
             color=ACCENT, bbox=dict(boxstyle="round", fc="white",
                                     ec="#cccccc", alpha=0.9))

    fig.suptitle("Reaction–diffusion coupling: mass-conserving rate enhancement "
                 "inside condensates", fontsize=14.5, y=1.02)
    fig.tight_layout()
    out = FIGURES / "nld_reaction_coupling.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"  saved {out}")


# ===========================================================================
# Hero composite — snapshots + key curves in one PRL-style figure
# ===========================================================================

def plot_hero(sim, rec, snaps, alpha):
    set_style()
    fig = plt.figure(figsize=(15, 8.6))
    gs = fig.add_gridspec(3, 4, height_ratios=[1, 1, 1.1],
                          hspace=0.42, wspace=0.32)

    # top row: coarsening sequence — pick frames nearest these target times
    targets = [20.0, 150.0, 600.0, 1500.0]
    show = [min(snaps, key=lambda s: abs(s[0] - tt)) for tt in targets]
    for i, (t, phi) in enumerate(show):
        ax = fig.add_subplot(gs[0, i])
        ax.imshow(phi, extent=[0, sim.L, 0, sim.L], cmap=DIV_CMAP,
                  vmin=-1, vmax=1, origin="lower", interpolation="bilinear")
        ax.set_title(fr"$t={t:.0f}$", fontsize=11)
        ax.set_xticks([]); ax.set_yticks([])

    # dispersion (mid-left) — Flory–Huggins thermodynamics
    axd = fig.add_subplot(gs[1, 0])
    k = np.linspace(0, 3, 400)
    chi_fh, N1_fh, N2_fh, kap_fh = 1.5, 100, 1, 0.5
    for phi0, c in [(0.10, ACCENT2), (0.30, ACCENT)]:
        fpp = d2f_dphi2(phi0, chi_fh, N1_fh, N2_fh)
        axd.plot(k, dispersion(k, fpp, 1.0, kap_fh), color=c, lw=2,
                 label=fr"$\phi_0={phi0}$")
        ks, ws = fastest_mode(fpp, kap_fh, 1.0)
        if np.isfinite(ks):
            axd.plot(ks, ws, 'o', color=c, ms=5)
    axd.axhline(0, color=INK, lw=0.7)
    axd.set_title(r"$\omega(k)$  (Flory–Huggins)", fontsize=11)
    axd.set_xlabel("$k$"); axd.legend(fontsize=8)

    # coarsening
    axc = fig.add_subplot(gs[1, 1])
    m = rec["t"] > 2
    axc.loglog(rec["t"][m], rec["L"][m], 'o', color=ACCENT, ms=3)
    tt = rec["t"][m]
    cref = rec["L"][m][len(tt)//2] / tt[len(tt)//2]**(1/3)
    axc.loglog(tt, cref * tt**(1/3), '--', color=INK, lw=1.8)
    axc.set_title(fr"$L\sim t^{{{alpha:.2f}}}$", fontsize=11)
    axc.set_xlabel("$t$")

    # energy
    axe = fig.add_subplot(gs[1, 2])
    axe.plot(rec["t"], rec["F"], color=INK, lw=2)
    axe.plot(rec["t"], rec["F_bulk"], color=ACCENT2, lw=1.5)
    axe.plot(rec["t"], rec["F_grad"], color=ACCENT, lw=1.5)
    axe.set_title(r"$F[\phi]\downarrow$", fontsize=11)
    axe.set_xlabel("$t$")

    # phase portrait
    axp = fig.add_subplot(gs[1, 3])
    axp.plot(rec["dense_frac"], rec["F"], color=ACCENT, lw=2)
    axp.plot(rec["dense_frac"][-1], rec["F"][-1], '*', color=INK, ms=12)
    axp.set_title("phase portrait", fontsize=11)
    axp.set_xlabel("$p(t)$"); axp.set_ylabel("$F$")

    # bottom row spanning: reaction coupling + bifurcation
    axr = fig.add_subplot(gs[2, :2])
    pp = np.linspace(1e-4, 1, 400)
    cmap = mpl.colormaps["viridis"]
    for i, Kv in enumerate([10, 100, 1000]):
        axr.plot(pp, effective_rate_multiplier(Kv, Kv, pp),
                 color=cmap(0.1 + 0.8 * i / 2), lw=2.2, label=fr"$K={Kv}$")
    axr.axhline(1, color=INK, ls=":", lw=0.8)
    axr.set_yscale("log")
    axr.set_xlabel(r"dense fraction $p$")
    axr.set_ylabel(r"$k_{\rm eff}/k_{\rm out}$")
    axr.set_title("Reaction enhancement (mass-conserving)", fontsize=11)
    axr.legend(fontsize=8)

    axb = fig.add_subplot(gs[2, 2:])
    N1, N2 = 100, 1
    phi_c, chi_c = critical_point(N1, N2)
    from flory_huggins import spinodal_curve
    phi_sp, chi_sp = spinodal_curve(N1, N2, 300)
    chi_bn, phL, phR = robust_binodal(N1, N2, chi_max_factor=3.0, n=150)
    if len(chi_bn):
        axb.plot(phL, chi_bn, color=ACCENT2, lw=2.2)
        axb.plot(phR, chi_bn, color=ACCENT2, lw=2.2)
    mid = len(phi_sp)//2
    axb.plot(phi_sp[:mid], chi_sp[:mid], '--', color=ACCENT, lw=1.8)
    axb.plot(phi_sp[mid:], chi_sp[mid:], '--', color=ACCENT, lw=1.8)
    axb.plot(phi_c, chi_c, 'o', color=INK, ms=6)
    axb.plot(0.3, 1.5, '*', color=GOLD, ms=14, mec=INK, mew=0.5)
    axb.set_xlim(0, 0.8); axb.set_ylim(0, chi_c*3)
    axb.set_xlabel(r"$\phi$"); axb.set_ylabel(r"$\chi$")
    axb.set_title("Phase diagram", fontsize=11)

    fig.suptitle("Nonlinear dynamics of biomolecular condensates: "
                 "Cahn–Hilliard spinodal decomposition + reaction coupling",
                 fontsize=16, y=0.995)
    out = FIGURES / "nld_hero_overview.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"  saved {out}")


# ===========================================================================
# Main
# ===========================================================================

def main():
    print("Nonlinear dynamics analysis of the Cahn-Hilliard + reaction model")
    print("=" * 66)

    # static (cheap) analyses
    print("\n[1] Linear stability / dispersion relation")
    plot_dispersion_relation()
    plot_fastest_mode_map()

    print("\n[2] Bifurcation diagram")
    plot_bifurcation_diagram()

    print("\n[6] Reaction-diffusion coupling")
    plot_reaction_coupling()

    # one expensive solver run feeds 3,4,5 and the hero
    print("\n[3-5] Cahn-Hilliard diagnostics (single run)")
    sim, rec, snaps = run_diagnostics()
    print("  Running Flory-Huggins solver for coarsening-arrest comparison...")
    fh_curve = fh_frozen_coarsening()
    alpha = plot_coarsening(rec, fh_curve=fh_curve)
    plot_energy_landscape(rec)
    plot_phase_portrait(rec, snaps)

    print("\n[hero] Composite overview")
    plot_hero(sim, rec, snaps, alpha)

    print("\nDone.")


if __name__ == "__main__":
    main()
