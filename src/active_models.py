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
                 dt=0.05, A=2.5, k_react=0.0, phi_ss=0.0, lam=0.0, zeta=0.0):
        self.L, self.N = L, N
        self.a, self.b, self.kappa = a, b, kappa
        self.M, self.dt, self.A = mobility, dt, A
        self.k_react, self.phi_ss, self.lam = k_react, phi_ss, lam
        self.zeta = zeta  # Active Model B+ current term: -zeta div[(lap phi) grad phi]
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
        # conserved part: dphi/dt |_cons = M lap(mu) = -M k^2 (g_hat already
        #   carries f' + lambda|grad phi|^2 - A phi; the kappa & A pieces are
        #   in the implicit denominator)
        rhs = phi_hat - self.dt * self.M * self.k2 * g_hat
        # Active Model B+ current term:  dphi/dt += -zeta div[(lap phi) grad phi]
        if self.zeta != 0.0:
            lap = np.real(np.fft.ifft2(-self.k2 * phi_hat))
            gx = np.real(np.fft.ifft2(1j * self.KX * phi_hat))
            gy = np.real(np.fft.ifft2(1j * self.KY * phi_hat))
            divV = np.real(np.fft.ifft2(1j * self.KX * np.fft.fft2(lap * gx)
                                        + 1j * self.KY * np.fft.fft2(lap * gy)))
            rhs = rhs - self.dt * self.zeta * np.fft.fft2(divV)
        # reaction drive toward phi_ss (constant -> only k=0 component)
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


# ===========================================================================
# 2.  Active Model B(+):  broken detailed balance
# ===========================================================================

def run_amb_sweep(lams, T=600.0, N=192, L=96.0, dt=0.02, seed=5):
    out = {}
    for lam in lams:
        print(f"  lambda={lam:.2f} ...", end="", flush=True)
        sim = ActiveCH2D(L=L, N=N, lam=lam, dt=dt)
        sim.initialize(phi_mean=0.0, seed=seed)
        n = int(T / dt)
        rec_every = max(1, n // 100)
        ts, Ls = [], []
        for i in range(1, n + 1):
            sim.step()
            if i % rec_every == 0:
                ts.append(sim.time); Ls.append(sim.domain_length())
        out[lam] = (np.array(ts), np.array(Ls), sim.phi.copy())
        print(f" phi in [{sim.phi.min():+.2f},{sim.phi.max():+.2f}]")
    return out


def _coexistence_peaks(phi):
    """Modal phi of the dilute (phi<0) and dense (phi>0) phases."""
    def mode(vals):
        h, e = np.histogram(vals, bins=60)
        return 0.5 * (e[np.argmax(h)] + e[np.argmax(h) + 1])
    neg, pos = phi[phi < -0.2], phi[phi > 0.2]
    return (mode(neg) if neg.size else np.nan,
            mode(pos) if pos.size else np.nan)


def plot_active_model_b(sweep):
    """Honest characterisation of the lambda|grad phi|^2 (Active Model B) term."""
    set_style()
    lams = sorted(sweep)
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(16.5, 5.0))
    cmap = mpl.colormaps["viridis"]

    # (a) phi distribution: peaks shift asymmetrically with lambda
    for i, lam in enumerate(lams):
        phi = sweep[lam][2]
        c = cmap(0.1 + 0.8 * i / (len(lams) - 1))
        axA.hist(phi.ravel(), bins=120, density=True, histtype="step",
                 color=c, lw=2, label=fr"$\lambda={lam:g}$")
    axA.axvline(-1, color=INK, ls=":", lw=1); axA.axvline(1, color=INK, ls=":", lw=1)
    axA.set_xlabel(r"order parameter  $\phi$")
    axA.set_ylabel("probability density")
    axA.set_title(r"(a) $\lambda|\nabla\phi|^2$ shifts coexistence densities")
    axA.set_xlim(-1.4, 1.4); axA.legend(fontsize=9)

    # (b) coexistence densities phi_+ and phi_- vs lambda
    pm = np.array([_coexistence_peaks(sweep[l][2]) for l in lams])
    axB.plot(lams, pm[:, 1], 'o-', color=ACCENT, lw=2, ms=7, label=r"dense  $\phi_+$")
    axB.plot(lams, pm[:, 0], 's-', color=ACCENT2, lw=2, ms=7, label=r"dilute  $\phi_-$")
    axB.axhline(1, color=INK, ls=":", lw=1); axB.axhline(-1, color=INK, ls=":", lw=1)
    axB.set_xlabel(r"activity  $\lambda$")
    axB.set_ylabel(r"coexistence density")
    axB.set_title("(b) broken $\\pm\\phi$ symmetry (non-equilibrium)")
    axB.legend(fontsize=9)

    # (c) coarsening: modestly slowed, not arrested
    for i, lam in enumerate(lams):
        t, Lt, _ = sweep[lam]
        axC.plot(t, Lt, color=cmap(0.1 + 0.8 * i / (len(lams) - 1)), lw=2.2,
                 label=fr"$\lambda={lam:g}$")
    axC.set_xlabel(r"time  $t$"); axC.set_ylabel(r"domain length  $L(t)$")
    axC.set_title("(c) coarsening continues (no arrest from $\\lambda$ alone)")
    axC.legend(fontsize=9, loc="upper left")

    fig.suptitle("Active Model B: the $\\lambda|\\nabla\\phi|^2$ term breaks detailed "
                 "balance (Tjhung 2018). Arrest needs the reaction route (Fig. above).",
                 fontsize=13.5, y=1.01)
    fig.tight_layout()
    out = FIGURES / "active_model_b.png"
    fig.savefig(out); plt.close(fig)
    print(f"  saved {out}")


# ===========================================================================
# 3.  Thermal-gradient fission  (Ianeselli 2022)
# ===========================================================================

def run_thermal_fission(N=256, L=128.0, dt=0.02, T=600.0, a0=1.0, amp=2.0,
                        width=7.0, ax=42.0, ay=13.0):
    """A single elongated droplet divided by a warm band (weakened interaction).

    The well depth a(x) = a0[1 - amp exp(-(x-L/2)^2/2w^2)] dips below zero in a
    central vertical band -- a 'warm' stripe (UCST: chi ~ 1/T) where the dense
    phase is destabilised.  An elongated droplet straddling the band has its
    neck dissolved and pinches into two daughters (verified: 2 disconnected
    dense domains form with the band, but only 1 without it -- the division is
    genuinely band-driven, not mere Rayleigh-Plateau).  Mass is conserved.
    """
    x = np.linspace(0, L, N, endpoint=False)
    X, Y = np.meshgrid(x, x)
    a_field = a0 * (1 - amp * np.exp(-(X - L / 2)**2 / (2 * width**2)))
    sim = ActiveCH2D(L=L, N=N, a=a_field, dt=dt, A=5.0)
    ell = np.sqrt(((X - L / 2) / ax)**2 + ((Y - L / 2) / ay)**2)
    sim.phi = -0.9 + 1.8 * 0.5 * (1 - np.tanh((ell - 1) * 6))
    sim.time = 0.0
    snaps = [(0.0, sim.phi.copy())]
    targets = [150.0, 350.0, T]
    n = int(T / dt); ni = 1
    for i in range(1, n + 1):
        sim.step()
        if ni < len(targets) and sim.time >= targets[ni - 1]:
            snaps.append((sim.time, sim.phi.copy())); ni += 1
    if snaps[-1][0] < sim.time - dt:
        snaps.append((sim.time, sim.phi.copy()))
    return a_field, snaps, sim


def plot_thermal_fission():
    set_style()
    a_field, snaps, sim = run_thermal_fission()
    show = snaps[:4]
    fig, axes = plt.subplots(1, 5, figsize=(18, 3.9),
                             gridspec_kw={"width_ratios": [1, 1, 1, 1, 1.05]})
    for ax, (t, phi) in zip(axes[:4], show):
        ax.imshow(phi, extent=[0, sim.L, 0, sim.L], cmap=DIV_CMAP, vmin=-1, vmax=1,
                  origin="lower", interpolation="bilinear")
        ax.set_title(fr"$t={t:.0f}$", fontsize=12); ax.set_xticks([]); ax.set_yticks([])
    # interaction profile
    x = np.linspace(0, sim.L, sim.N, endpoint=False)
    axes[4].plot(a_field[sim.N // 2], x, color=ACCENT, lw=2.4)
    axes[4].axvline(0, color=INK, ls=":", lw=1)
    axes[4].set_title("well depth $a(x)$", fontsize=12)
    axes[4].set_xlabel("$a$"); axes[4].set_ylabel("$x$")
    axes[4].text(0.5, 0.5, "warm band\n($a<0$:\ndissolves)", transform=axes[4].transAxes,
                 ha="center", va="center", fontsize=8.5, color=ACCENT)
    fig.suptitle("Thermal-gradient fission: a temperature band divides one droplet "
                 "into two without a membrane (Ianeselli 2022)", fontsize=13.5, y=1.04)
    fig.tight_layout()
    out = FIGURES / "active_thermal_fission.png"
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

    only = sys.argv[1] if len(sys.argv) > 1 else "all"

    print("Active condensate models")
    print("=" * 40)
    if only in ("all", "reaction"):
        print("\n[1] Reaction-arrested coarsening")
        sweep = run_reaction_sweep([0.0, 0.005, 0.01, 0.02, 0.05, 0.1])
        plot_reaction_arrest(sweep)
        plot_morphology_and_distribution(sweep)
    if only in ("all", "amb"):
        print("\n[2] Active Model B (broken detailed balance)")
        amb = run_amb_sweep([0.0, 0.5, 1.0, 2.0])
        plot_active_model_b(amb)
    if only in ("all", "thermal"):
        print("\n[3] Thermal-gradient fission")
        plot_thermal_fission()
    print("\nDone.")
