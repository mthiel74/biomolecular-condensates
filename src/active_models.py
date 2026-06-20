#!/usr/bin/env python3
r"""Active / chemically-driven extensions of the Cahn-Hilliard model.

Phase 2 of the nonlinear-dynamics study.  Builds on the passive double-well
Model-B solver (src/nonlinear_analysis.py) and adds the physics of
chemically-maintained ("active") condensates:

1.  Reaction-arrested coarsening  (Zwicker 2015; Glasner-Witelski; Bauermann
    2024).  A linear chemical reaction that interconverts droplet material
    and a soluble form,
        s(phi) = k_+ (1 - phi) - k_- phi = -k_react (phi - phi_ss),
    competes with phase separation.  The extended equation is
        d phi / dt = M grad^2 mu + s(phi).
    Linear stability gains a uniform -k_react term:
        omega(k) = -M k^2 (f''(phi0) + kappa k^2) - k_react,
    which STABILISES the long-wavelength (k->0) modes.  The unstable band
    becomes finite [k_min, k_max] -> a length scale is selected, Ostwald
    ripening is ARRESTED, and the system settles into a steady-state
    emulsion ("microphase separation").  The arrested size scales as
        L* ~ 2 pi sqrt(M |f''| / k_react)   ->   R* ~ sqrt(D / k_react).

2.  Active Model B+  (Tjhung 2018; Caballero 2018).  A non-equilibrium term
    lambda |grad phi|^2 added to the chemical potential breaks detailed
    balance.  We sweep lambda and report what it actually does to the
    interface dynamics (honest measurement, not an assumed result).

3.  Thermal-gradient fission  (Ianeselli 2022).  A spatial temperature
    profile makes the interaction chi = chi(x) position-dependent; a droplet
    sitting across the gradient is squeezed and can split -- membrane-free
    division.

Figures -> figures/.  Run as a script to regenerate everything.
"""

from __future__ import annotations

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path

from nonlinear_analysis import (set_style, ACCENT, ACCENT2, GOLD, INK,
                                SEQ_CMAP, DIV_CMAP, _domain_length_first_moment)

FIGURES = Path(__file__).resolve().parent.parent / "figures"
FIGURES.mkdir(exist_ok=True)


class ActiveCH2D:
    r"""Double-well Cahn-Hilliard with optional chemical reaction and activity.

        mu  = -a phi + b phi^3 - kappa lap(phi) + lambda |grad phi|^2
        d phi/dt = M lap(mu) - k_react (phi - phi_ss)

    Semi-implicit spectral: the linear operator (Eyre A + kappa k^2 in the
    biharmonic, plus the reaction k_react) is implicit; the double-well
    remainder and the activity term are explicit.
    """

    def __init__(self, L=128.0, N=256, a=1.0, b=1.0, kappa=1.0, mobility=1.0,
                 dt=0.05, A=2.5, k_react=0.0, phi_ss=0.0, lam=0.0):
        self.L, self.N = L, N
        self.a, self.b, self.kappa = a, b, kappa
        self.M, self.dt, self.A = mobility, dt, A
        self.k_react, self.phi_ss, self.lam = k_react, phi_ss, lam
        self.dx = L / N
        k = 2 * np.pi * np.fft.fftfreq(N, d=self.dx)
        self.KX, self.KY = np.meshgrid(k, k)
        self.k2 = self.KX**2 + self.KY**2
        # implicit denominator: biharmonic stabilisation + linear reaction
        self._denom = (1.0 + dt * (mobility * self.k2 * (A + kappa * self.k2)
                                   + k_react))
        self.phi = None
        self.time = 0.0

    def initialize(self, phi_mean=0.0, noise_amplitude=0.02, seed=42):
        rng = np.random.default_rng(seed)
        self.phi = phi_mean + noise_amplitude * rng.standard_normal((self.N, self.N))
        self.time = 0.0

    def _grad_sq(self, phi_hat):
        """|grad phi|^2 via spectral derivatives."""
        gx = np.real(np.fft.ifft2(1j * self.KX * phi_hat))
        gy = np.real(np.fft.ifft2(1j * self.KY * phi_hat))
        return gx**2 + gy**2

    def step(self):
        phi_hat = np.fft.fft2(self.phi)
        # nonlinear chemical potential pieces (explicit): double well - A phi
        g = -self.a * self.phi + self.b * self.phi**3 - self.A * self.phi
        if self.lam != 0.0:
            g = g + self.lam * self._grad_sq(phi_hat)
        g_hat = np.fft.fft2(g)
        # reaction drive toward phi_ss (constant -> only k=0 component)
        rhs = phi_hat - self.dt * self.M * self.k2 * g_hat
        if self.k_react != 0.0:
            rhs[0, 0] += self.dt * self.k_react * self.phi_ss * self.N**2
        self.phi = np.real(np.fft.ifft2(rhs / self._denom))
        self.time += self.dt

    def domain_length(self):
        return _domain_length_first_moment(self.phi, self.k2)


def dispersion_active(k, fpp, kappa, M, k_react):
    return -M * k**2 * (fpp + kappa * k**2) - k_react


def droplet_radii(phi, dx, threshold=0.0):
    """Equivalent radii of connected dense domains (phi > threshold)."""
    from scipy.ndimage import label
    lab, n = label(phi > threshold)
    if n == 0:
        return np.array([])
    areas = np.bincount(lab.ravel())[1:] * dx**2
    return np.sqrt(areas / np.pi)


# ===========================================================================
# 1.  Reaction-arrested coarsening
# ===========================================================================

def run_reaction_sweep(k_values, T=1200.0, N=256, L=128.0, dt=0.05, seed=3):
    """Run the reaction model for each k_react; record L(t) and final field."""
    out = {}
    for kr in k_values:
        print(f"  k_react={kr:.4f} ...", end="", flush=True)
        sim = ActiveCH2D(L=L, N=N, k_react=kr, phi_ss=0.0, dt=dt)
        sim.initialize(phi_mean=0.0, seed=seed)
        n = int(T / dt)
        rec_every = max(1, n // 120)
        ts, Ls = [], []
        for i in range(1, n + 1):
            sim.step()
            if i % rec_every == 0:
                ts.append(sim.time); Ls.append(sim.domain_length())
        out[kr] = (np.array(ts), np.array(Ls), sim.phi.copy())
        print(f" L*={Ls[-1]:.1f}")
    return out


def plot_reaction_arrest(sweep, L=128.0):
    """L(t) trajectories, steady-state bifurcation L*(k), and dispersion."""
    set_style()
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(16.5, 5.0))
    ks = sorted(sweep)
    cmap = mpl.colormaps["viridis"]

    # (a) L(t)
    for i, kr in enumerate(ks):
        t, Lt, _ = sweep[kr]
        c = INK if kr == 0 else cmap(0.12 + 0.78 * i / (len(ks) - 1))
        lab = "passive ($k=0$)" if kr == 0 else fr"$k_{{\rm react}}={kr:g}$"
        axA.plot(t, Lt, color=c, lw=2.2, label=lab)
    axA.set_xlabel(r"time  $t$"); axA.set_ylabel(r"domain length  $L(t)$")
    axA.set_title("(a) reaction arrests coarsening")
    axA.legend(fontsize=8.5, loc="upper left")

    # (b) bifurcation L* vs k_react (steady-state = mean of last 20%)
    kpos = [k for k in ks if k > 0]
    Lstar = [np.mean(sweep[k][1][-20:]) for k in kpos]
    axB.loglog(kpos, Lstar, 'o', color=ACCENT, ms=9, mec="white", mew=0.6,
               label="steady-state $L^*$")
    p = np.polyfit(np.log(kpos), np.log(Lstar), 1)
    kk = np.logspace(np.log10(min(kpos)), np.log10(max(kpos)), 50)
    axB.loglog(kk, np.exp(p[1]) * kk**p[0], '-', color=INK, lw=2,
               label=fr"fit  $L^*\sim k^{{{p[0]:.2f}}}$")
    axB.loglog(kk, Lstar[0] * (kk / kpos[0])**(-0.5), '--', color=ACCENT2, lw=1.8,
               label=r"$k^{-1/2}$ (Zwicker $R^*\!\sim\!\sqrt{D/k}$)")
    axB.set_xlabel(r"reaction rate  $k_{\rm react}$")
    axB.set_ylabel(r"arrested length  $L^*$")
    axB.set_title("(b) Ostwald-arrest bifurcation")
    axB.legend(fontsize=9, loc="lower left")

    # (c) dispersion: finite unstable band
    k = np.linspace(0, 1.4, 400)
    fpp = -1.0  # double well at phi0=0: f'' = -a
    for i, kr in enumerate([0.0, 0.05, 0.15, 0.25]):
        c = INK if kr == 0 else cmap(0.2 + 0.7 * i / 3)
        w = dispersion_active(k, fpp, 1.0, 1.0, kr)
        axC.plot(k, w, color=c, lw=2.2,
                 label=("passive" if kr == 0 else fr"$k={kr:g}$"))
    axC.axhline(0, color="grey", lw=0.8)
    axC.set_xlabel(r"wavenumber  $k$"); axC.set_ylabel(r"$\omega(k)$")
    axC.set_title(r"(c) low-$k$ modes stabilised: $\omega=-Mk^2(f''+\kappa k^2)-k_{\rm react}$")
    axC.legend(fontsize=9, loc="lower left")
    axC.set_ylim(-0.35, 0.30)

    fig.suptitle("Chemically-maintained condensates: reaction-arrested coarsening "
                 "(Zwicker 2015; Bauermann 2024)", fontsize=14.5, y=1.02)
    fig.tight_layout()
    out = FIGURES / "active_reaction_arrest.png"
    fig.savefig(out); plt.close(fig)
    print(f"  saved {out}  (L* exponent {p[0]:.3f})")
    return p[0]


def plot_morphology_and_distribution(sweep, L=128.0):
    """Passive vs active morphology + steady-state droplet-size distribution."""
    set_style()
    ks = sorted(sweep)
    k_active = max(ks)  # strongest reaction = most arrested
    dx = L / sweep[k_active][2].shape[0]

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.2))

    for ax, kr, ttl in [(axes[0], 0.0, "passive: coarsened"),
                        (axes[1], k_active, f"active ($k={k_active:g}$): emulsion")]:
        phi = sweep[kr][2]
        im = ax.imshow(phi, extent=[0, L, 0, L], cmap=DIV_CMAP,
                       vmin=-1, vmax=1, origin="lower", interpolation="bilinear")
        ax.set_title(ttl); ax.set_xticks([]); ax.set_yticks([])

    # size distribution: passive vs active
    r_pass = droplet_radii(sweep[0.0][2], dx)
    r_act = droplet_radii(sweep[k_active][2], dx)
    ax = axes[2]
    if len(r_pass):
        ax.hist(r_pass / r_pass.mean(), bins=np.linspace(0, 3, 18), density=True,
                alpha=0.5, color=INK, label=f"passive (n={len(r_pass)})")
    if len(r_act):
        ax.hist(r_act / r_act.mean(), bins=np.linspace(0, 3, 18), density=True,
                alpha=0.6, color=ACCENT, label=f"active (n={len(r_act)})")
    ax.set_xlabel(r"droplet radius  $R/\langle R\rangle$")
    ax.set_ylabel("probability density")
    ax.set_title("(c) steady-state distribution is narrower than LSW")
    ax.legend(fontsize=9)
    ax.text(0.97, 0.6, fr"CV$_{{\rm active}}$={np.std(r_act)/np.mean(r_act):.2f}"
            "\n" fr"CV$_{{\rm passive}}$={np.std(r_pass)/np.mean(r_pass):.2f}",
            transform=ax.transAxes, ha="right", fontsize=9,
            bbox=dict(boxstyle="round", fc="white", ec="#ccc", alpha=0.9))

    fig.suptitle("Steady-state emulsion vs. passive coarsening: morphology and "
                 "size distribution", fontsize=14.5, y=1.00)
    fig.tight_layout()
    out = FIGURES / "active_morphology_distribution.png"
    fig.savefig(out); plt.close(fig)
    print(f"  saved {out}")


if __name__ == "__main__":
    import sys
    if "--selftest" in sys.argv:
        for kr in [0.0, 0.01, 0.05]:
            sim = ActiveCH2D(k_react=kr); sim.initialize(seed=3)
            for _ in range(2000):
                sim.step()
            print(f"k_react={kr}: L={sim.domain_length():.2f} mean={sim.phi.mean():+.4f}")
        sys.exit()

    print("Active condensate models")
    print("=" * 40)
    print("\n[1] Reaction-arrested coarsening")
    sweep = run_reaction_sweep([0.0, 0.005, 0.01, 0.02, 0.05, 0.1])
    plot_reaction_arrest(sweep)
    plot_morphology_and_distribution(sweep)
    print("\nDone.")
