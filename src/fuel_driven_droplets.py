#!/usr/bin/env python3
r"""Fuel-driven droplets and the arrest of Ostwald ripening
(Zwicker, Hyman & Jülicher 2015, Nat. Phys.; Weber, Zwicker, Jülicher & Lee
2019, Rep. Prog. Phys.; Lee 2021).

Phase 2 of the nonlinear-dynamics study.  A chemical reaction that
interconverts droplet material (dense phase) and a soluble form (dilute
phase) is added as a source term to the conserved Cahn-Hilliard dynamics:

    d phi/dt = M lap(mu) + s(phi),    mu = f'(phi) - kappa lap(phi)
    s(phi)   = k_+ (1 - phi) - k_- (phi + 1)  =  -k_react (phi - phi_ss)

with k_react = k_+ + k_-  the total turnover rate and
phi_ss = (k_+ - k_-)/(k_+ + k_-) the composition the reaction drives toward.
(The double-well solver uses the symmetric convention phi in [-1, +1] with
dilute = -1, dense = +1; the user's [0,1] form maps onto this linear
relaxation toward phi_ss.)

The reaction continuously *pumps* material between the phases, holding the
system away from equilibrium and qualitatively changing the physics:

  * Linear stability gains a uniform -k_react term:
        omega(k) = -M k^2 (f''(phi0) + kappa k^2) - k_react.
    The long-wavelength (k -> 0) modes that drive Ostwald ripening are now
    STABILISED.  The unstable band becomes finite, [k_min, k_max], so a
    length scale is selected and ripening is ARRESTED.

  * There is a hard DISSOLUTION threshold: the fastest passive mode grows at
    omega_max = M f''^2 / (4 kappa); once k_react exceeds it,
        k_diss = M f''(phi0)^2 / (4 kappa),
    every mode is stable and the homogeneous (mixed) state wins -- droplets
    dissolve.  So along k_react one passes
        Ostwald (k->0)  ->  arrested emulsion (0 < k < k_diss)  ->  dissolution.

  * The arrested size scales as  L* ~ 2 pi sqrt(M |f''| / k_react)  ->
    R* ~ sqrt(D / k_react)  (Zwicker 2015), and many droplets coexist at
    steady state instead of one big droplet winning.

Module E (Lee 2021): in the conversion-limited transient the coarsening
exponent is measured and compared with the passive Lifshitz-Slyozov t^{1/3}.

Figures -> figures/.  Run as a script to regenerate everything:
    python3 src/fuel_driven_droplets.py
"""

from __future__ import annotations

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path

from nonlinear_analysis import (set_style, ACCENT, ACCENT2, GOLD, INK,
                                SEQ_CMAP, DIV_CMAP)

FIGURES = Path(__file__).resolve().parent.parent / "figures"
FIGURES.mkdir(exist_ok=True)


# ===========================================================================
# Solver
# ===========================================================================

class FuelDrivenCH2D:
    r"""Double-well Cahn-Hilliard with a linear chemical reaction source.

        mu = -a phi + b phi^3 - kappa lap(phi)
        d phi/dt = M lap(mu) - k_react (phi - phi_ss)

    Semi-implicit spectral: the Eyre-stabilised biharmonic AND the linear
    reaction are treated implicitly (the reaction adds +dt*k_react to the
    denominator, which is unconditionally stabilising); the double-well
    remainder is explicit.
    """

    def __init__(self, L=128.0, N=256, a=1.0, b=1.0, kappa=1.0, mobility=1.0,
                 dt=0.05, A=2.5, k_react=0.0, phi_ss=0.0):
        self.L, self.N = L, N
        self.a, self.b, self.kappa = a, b, kappa
        self.M, self.dt, self.A = mobility, dt, A
        self.k_react, self.phi_ss = k_react, phi_ss
        self.dx = L / N
        k = 2 * np.pi * np.fft.fftfreq(N, d=self.dx)
        self.KX, self.KY = np.meshgrid(k, k)
        self.k2 = self.KX**2 + self.KY**2
        self._denom = (1.0 + dt * (mobility * self.k2 * (A + kappa * self.k2)
                                   + k_react))
        self.phi = None
        self.time = 0.0

    def initialize(self, phi_mean=0.0, noise_amplitude=0.05, seed=42):
        rng = np.random.default_rng(seed)
        self.phi = phi_mean + noise_amplitude * rng.standard_normal((self.N, self.N))
        self.time = 0.0

    def step(self):
        phi_hat = np.fft.fft2(self.phi)
        g = (-self.a * self.phi + self.b * self.phi**3) - self.A * self.phi
        g_hat = np.fft.fft2(g)
        rhs = phi_hat - self.dt * self.M * self.k2 * g_hat
        # reaction drive toward phi_ss enters only the k=0 (mean) component
        if self.k_react != 0.0:
            rhs[0, 0] += self.dt * self.k_react * self.phi_ss * self.N**2
        self.phi = np.real(np.fft.ifft2(rhs / self._denom))
        self.time += self.dt

    def run(self, n_steps):
        for _ in range(n_steps):
            self.step()

    def domain_length(self):
        """L = 2 pi <S> / <k S> from the structure factor (no peak-finding)."""
        phi = self.phi - self.phi.mean()
        S = np.abs(np.fft.fft2(phi))**2
        kmag = np.sqrt(self.k2)
        den = np.sum(kmag * S)
        return 2 * np.pi * np.sum(S) / den if den > 0 else self.L

    def order_param(self):
        """RMS order parameter std(phi): -> 0 when droplets dissolve."""
        return float(self.phi.std())

    def dense_fraction(self, threshold=0.5):
        """Fraction of pixels in the well-separated dense phase (phi > 0.5)."""
        return float(np.mean(self.phi > threshold))


def dispersion(k, fpp, kappa, M, k_react):
    return -M * k**2 * (fpp + kappa * k**2) - k_react


def k_dissolution(fpp, kappa, M):
    """Reaction rate above which every mode is stable (droplets dissolve)."""
    return M * fpp**2 / (4.0 * kappa)


def band_edges(k_react, fpp, kappa, M):
    """Exact edges of the unstable band, roots of omega(k)=0.

    Mkappa k^4 - M|fpp| k^2 + k_react = 0  ->  two k^2 roots.
    Returns (k_min, k_max) or (nan, nan) above the dissolution threshold.
    """
    a = -fpp                                   # |fpp| > 0 inside the spinodal
    disc = (M * a)**2 - 4 * M * kappa * k_react
    if disc < 0:
        return np.nan, np.nan
    s = np.sqrt(disc)
    k2_min = (M * a - s) / (2 * M * kappa)
    k2_max = (M * a + s) / (2 * M * kappa)
    return np.sqrt(max(k2_min, 0.0)), np.sqrt(max(k2_max, 0.0))


def fastest_mode_wavelength(fpp, kappa):
    """lambda* = 2 pi / k*, k* = sqrt(-fpp/2kappa): the k-independent floor."""
    return 2 * np.pi / np.sqrt(-fpp / (2 * kappa))


def droplet_radii(phi, dx, threshold=0.0):
    """Equivalent radii of connected dense domains (phi > threshold)."""
    from scipy.ndimage import label
    lab, n = label(phi > threshold)
    if n == 0:
        return np.array([])
    areas = np.bincount(lab.ravel())[1:] * dx**2
    return np.sqrt(areas / np.pi)


# ===========================================================================
# Sweeps
# ===========================================================================

def run_rate_sweep(k_values, T=1400.0, N=256, L=128.0, dt=0.05,
                   phi_mean=0.0, phi_ss=0.0, seed=3, n_rec=130):
    """Evolve the model for each k_react; record L(t), field, dense fraction."""
    out = {}
    for kr in k_values:
        print(f"  k_react={kr:.4f} ...", end="", flush=True)
        sim = FuelDrivenCH2D(L=L, N=N, k_react=kr, phi_ss=phi_ss, dt=dt)
        sim.initialize(phi_mean=phi_mean, seed=seed)
        n = int(T / dt)
        rec_every = max(1, n // n_rec)
        ts, Ls = [], []
        for i in range(1, n + 1):
            sim.step()
            if i % rec_every == 0:
                ts.append(sim.time)
                Ls.append(sim.domain_length())
        out[kr] = (np.array(ts), np.array(Ls), sim.phi.copy(),
                   sim.order_param())
        print(f" L*={Ls[-1]:.1f}  std(phi)={out[kr][3]:.3f}", flush=True)
    return out


# ===========================================================================
# Figures
# ===========================================================================

def plot_bifurcation(sweep, fpp=-1.0, kappa=1.0, M=1.0,
                     fname="fuel_bifurcation.png"):
    """Three-panel: L*(k) bracketed by analytic band edges, dissolution of the
    order parameter, and the dispersion relation with its finite band."""
    set_style()
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(16.5, 5.0))
    ks = sorted(sweep)
    kpos = np.array([k for k in ks if k > 0])
    kdiss = k_dissolution(fpp, kappa, M)
    lam_star = fastest_mode_wavelength(fpp, kappa)

    # (a) bifurcation L*(k): data, fit, k^-1/2 guide, floor, band-edge envelope
    Lstar = np.array([np.mean(sweep[k][1][-12:]) for k in kpos])
    op = np.array([sweep[k][3] for k in kpos])
    alive = op > 0.05                      # genuine phase separation present
    kk = np.logspace(np.log10(min(kpos) * 0.8), np.log10(kdiss * 0.999), 200)
    lam_max = np.array([2 * np.pi / band_edges(k, fpp, kappa, M)[0] for k in kk])
    # measured data
    axA.loglog(kpos[alive], Lstar[alive], 'o', color=ACCENT, ms=9,
               mec="white", mew=0.6, zorder=6, label=r"steady-state $L^*$ (2D sim)")
    # power-law fit through the data
    p = np.polyfit(np.log(kpos[alive]), np.log(Lstar[alive]), 1)
    axA.loglog(kk, np.exp(p[1]) * kk**p[0], '-', color=INK, lw=1.8,
               label=fr"fit $L^*\sim k^{{{p[0]:.2f}}}$")
    # k^-1/2 guide anchored to the small-k datum (Zwicker R*~sqrt(D/k))
    k0, L0 = kpos[alive][0], Lstar[alive][0]
    axA.loglog(kk, L0 * (kk / k0)**(-0.5), ':', color=GOLD, lw=1.9,
               label=r"$k^{-1/2}$ (Zwicker $R^*\!\sim\!\sqrt{D/k}$)")
    # theoretical envelope: floor lambda* and band edge lambda_max
    axA.loglog(kk, np.full_like(kk, lam_star), '--', color=ACCENT2, lw=1.4,
               label=r"$\lambda^*$ floor (fastest mode)")
    axA.loglog(kk, lam_max, '-.', color=ACCENT2, lw=1.2, alpha=0.7,
               label=r"$\lambda_{\max}$ band edge")
    axA.axvline(kdiss, color="grey", lw=1.4, ls=":")
    axA.text(kdiss * 0.96, lam_star * 1.04, r"$k_{\rm diss}$", color="grey",
             fontsize=10, ha="right", va="bottom", rotation=90)
    axA.set_xlabel(r"reaction rate  $k_{\rm react}$")
    axA.set_ylabel(r"arrested length scale")
    axA.set_title("(a) reaction selects a length (Ostwald arrest)")
    axA.legend(fontsize=7.5, loc="lower left")

    # (b) order parameter std(phi) vs k: collapses to zero at dissolution
    axB.plot(kpos, op, 'o-', color=ACCENT, ms=8, mec="white", mew=0.6, lw=2)
    axB.axvline(kdiss, color=GOLD, lw=1.8, ls=":")
    axB.fill_betweenx([0, 1], kdiss, max(kpos) * 1.5, color=GOLD, alpha=0.10)
    axB.text(kdiss * 1.04, op.max() * 0.55, "dissolution\n(mixed state)",
             color="#8a6a1a", fontsize=10, ha="left", va="center")
    axB.text(min(kpos) * 1.1, op.max() * 0.45, "arrested\nemulsion",
             color=ACCENT, fontsize=10, ha="left", va="center")
    axB.set_xlabel(r"reaction rate  $k_{\rm react}$")
    axB.set_ylabel(r"order parameter  $\langle\phi^2\rangle^{1/2}$")
    axB.set_title(r"(b) droplets dissolve above $k_{\rm diss}=M f''^2/4\kappa$")
    axB.set_ylim(0, op.max() * 1.15)
    axB.set_xlim(0, max(kpos) * 1.25)

    # (c) dispersion: finite unstable band + dissolution
    k = np.linspace(0, 1.4, 400)
    for i, kr in enumerate([0.0, 0.05, 0.15, kdiss]):
        c = INK if kr == 0 else mpl.colormaps["viridis"](0.2 + 0.7 * i / 3)
        w = dispersion(k, fpp, kappa, M, kr)
        lab = ("passive" if kr == 0 else
               (r"$k_{\rm diss}$" if abs(kr - kdiss) < 1e-9 else fr"$k={kr:g}$"))
        axC.plot(k, w, color=c, lw=2.2, label=lab)
    axC.axhline(0, color="grey", lw=0.8)
    axC.set_xlabel(r"wavenumber  $k$"); axC.set_ylabel(r"$\omega(k)$")
    axC.set_title(r"(c) $\omega=-Mk^2(f''+\kappa k^2)-k_{\rm react}$")
    axC.legend(fontsize=9, loc="lower left")
    axC.set_ylim(-0.32, 0.30)

    fig.suptitle("Fuel-driven droplets: chemical reactions arrest Ostwald "
                 "ripening  (Zwicker 2015; Weber 2019)", fontsize=14.5, y=1.02)
    fig.tight_layout()
    out = FIGURES / fname
    fig.savefig(out)
    plt.close(fig)
    print(f"  saved {out}")


def plot_passive_vs_fuel(sweep_droplet, L=128.0,
                         fname="fuel_passive_vs_active.png"):
    """Side-by-side morphology (passive Ostwald vs fuel-driven emulsion) and
    the steady-state droplet-size distribution vs the LSW form."""
    set_style()
    ks = sorted(sweep_droplet)
    k_passive = 0.0
    k_active = max(ks)
    dx = L / sweep_droplet[k_active][2].shape[0]

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.2))
    for ax, kr, ttl in [
        (axes[0], k_passive, "(a) passive  $k=0$:\none droplet wins (Ostwald)"),
        (axes[1], k_active, f"(b) fuel-driven  $k={k_active:g}$:\n"
                            "many droplets coexist")]:
        phi = sweep_droplet[kr][2]
        ax.imshow(phi, extent=[0, L, 0, L], cmap=DIV_CMAP, vmin=-1.1, vmax=1.1,
                  origin="lower", interpolation="bilinear")
        ax.set_title(ttl, fontsize=12)
        ax.set_xticks([]); ax.set_yticks([])

    # size distribution vs LSW
    ax = axes[2]
    r_act = droplet_radii(sweep_droplet[k_active][2], dx)
    if len(r_act):
        u = r_act / r_act.mean()
        ax.hist(u, bins=np.linspace(0, 2.5, 16), density=True, alpha=0.7,
                color=ACCENT, label=f"fuel-driven (n={len(r_act)})")
        cv = np.std(r_act) / np.mean(r_act)
        ax.text(0.96, 0.78, fr"CV $=$ {cv:.2f}", transform=ax.transAxes,
                ha="right", fontsize=10,
                bbox=dict(boxstyle="round", fc="white", ec="#ccc", alpha=0.9))
    # LSW 2D-style distribution (Lifshitz-Slyozov), broad with sharp cutoff at u=1.5
    u = np.linspace(0, 1.5, 300)
    g = np.where(u < 1.5,
                 (4 / 9.) * u**2 * (3 / (3 + u))**(7 / 3.)
                 * (1.5 / (1.5 - u))**(11 / 3.) * np.exp(-u / (1.5 - u)),
                 0.0)
    g = np.nan_to_num(g)
    if g.max() > 0:
        ax.plot(u, g / np.trapezoid(g, u), '-', color=INK, lw=2.4,
                label="LSW (Ostwald) distribution")
    ax.set_xlabel(r"droplet radius  $R/\langle R\rangle$")
    ax.set_ylabel("probability density")
    ax.set_title("(c) fuel-driven distribution is narrower than LSW")
    ax.legend(fontsize=9, loc="upper left")
    ax.set_xlim(0, 2.5)

    fig.suptitle("Passive coarsening vs. fuel-driven steady-state emulsion",
                 fontsize=14.5, y=1.00)
    fig.tight_layout()
    out = FIGURES / fname
    fig.savefig(out)
    plt.close(fig)
    print(f"  saved {out}")


# ===========================================================================
# Module E: coarsening exponents -- diffusion-limited vs conversion-limited
#
# Honest physics.  The fuel-driven model above uses a *bulk* reaction source,
# which ARRESTS coarsening (Module B) rather than producing a sustained power
# law.  The modified exponent of Lee (2021) belongs to the CONVERSION-LIMITED
# regime, where growth is controlled by the interface conversion reaction, not
# by long-range diffusion.  That regime is in the t^{1/2} universality class
# (Wagner 1961; Bray review) -- the same class as non-conserved curvature-
# driven (Allen-Cahn) coarsening -- whereas the passive, diffusion-limited
# conserved dynamics (Model B) gives the Lifshitz-Slyozov t^{1/3}.  We measure
# BOTH exponents directly and contrast them.
# ===========================================================================

class AllenCahn2D:
    r"""Non-conserved (Model A) double-well solver: the conversion/interface-
    limited universality class, L ~ t^{1/2}.

        d phi/dt = -M ( f'(phi) - kappa lap(phi) ),   f'(phi) = -a phi + b phi^3
    """

    def __init__(self, L=256.0, N=256, a=1.0, b=1.0, kappa=1.0, mobility=1.0,
                 dt=0.05, A=2.5):
        self.L, self.N, self.a, self.b, self.kappa = L, N, a, b, kappa
        self.M, self.dt, self.A = mobility, dt, A
        self.dx = L / N
        k = 2 * np.pi * np.fft.fftfreq(N, d=self.dx)
        KX, KY = np.meshgrid(k, k)
        self.k2 = KX**2 + KY**2
        self._denom = 1.0 + dt * mobility * (A + kappa * self.k2)
        self.phi = None
        self.time = 0.0

    def initialize(self, phi_mean=0.0, noise_amplitude=0.1, seed=11):
        rng = np.random.default_rng(seed)
        self.phi = phi_mean + noise_amplitude * rng.standard_normal((self.N, self.N))
        self.time = 0.0

    def step(self):
        g = (-self.a * self.phi + self.b * self.phi**3) - self.A * self.phi
        phi_hat = np.fft.fft2(self.phi)
        self.phi = np.real(np.fft.ifft2(
            (phi_hat - self.dt * self.M * np.fft.fft2(g)) / self._denom))
        self.time += self.dt

    def domain_length(self):
        phi = self.phi - self.phi.mean()
        S = np.abs(np.fft.fft2(phi))**2
        den = np.sum(np.sqrt(self.k2) * S)
        return 2 * np.pi * np.sum(S) / den if den > 0 else self.L


def plot_coarsening_exponent(fname="fuel_coarsening_exponent.png",
                             N=256, L=256.0, dt=0.05):
    """Diffusion-limited (conserved, t^1/3) vs conversion-limited
    (non-conserved, t^1/2) coarsening, both measured and fit."""
    from nonlinear_analysis import DoubleWellCH2D
    set_style()
    fig, ax = plt.subplots(figsize=(8.4, 6.2))

    def Lof(sim):
        phi = sim.phi - sim.phi.mean()
        S = np.abs(np.fft.fft2(phi))**2
        den = np.sum(np.sqrt(sim.k2) * S)
        return 2 * np.pi * np.sum(S) / den if den > 0 else sim.L

    def run(sim, T):
        n = int(T / dt)
        rec = max(1, n // 120)
        ts, Ls = [], []
        for i in range(1, n + 1):
            sim.step()
            if i % rec == 0 and sim.time > 2.0:
                ts.append(sim.time); Ls.append(Lof(sim))
        return np.array(ts), np.array(Ls)

    def fit_tail(t, Lc, frac=0.4):
        m = t > t[int((1 - frac) * len(t))]
        return np.polyfit(np.log(t[m]), np.log(Lc[m]), 1)[0]

    # conserved, diffusion-limited (Model B) -> t^1/3
    simB = DoubleWellCH2D(L=L, N=N, dt=dt, A=2.5)
    simB.initialize(phi_mean=0.0, noise_amplitude=0.05, seed=11)
    tB, LB = run(simB, 700.0)
    eB = fit_tail(tB, LB)
    ax.loglog(tB, LB, 'o', color=ACCENT2, ms=4.5,
              label="conserved, diffusion-limited (Model B)")
    tg = np.logspace(np.log10(tB[len(tB)//4]), np.log10(tB[-1]), 40)
    cB = LB[len(LB)//2] / tB[len(tB)//2]**(1/3)
    ax.loglog(tg, cB * tg**(1/3), '-', color=ACCENT2, lw=1.7,
              label=fr"fit $L\sim t^{{{eB:.2f}}}$  (LSW $1/3$)")

    # non-conserved, conversion/interface-limited (Model A) -> t^1/2
    simA = AllenCahn2D(L=L, N=N, dt=dt, A=2.5)
    simA.initialize(phi_mean=0.0, noise_amplitude=0.1, seed=11)
    tA, LA = run(simA, 700.0)
    eA = fit_tail(tA, LA)
    ax.loglog(tA, LA, 's', color=ACCENT, ms=4.5,
              label="conversion/interface-limited (Model A)")
    cA = LA[len(LA)//2] / tA[len(tA)//2]**0.5
    ax.loglog(tg, cA * tg**0.5, '-', color=ACCENT, lw=1.7,
              label=fr"fit $L\sim t^{{{eA:.2f}}}$  (Wagner/Lee $1/2$)")

    ax.set_xlabel(r"time  $t$")
    ax.set_ylabel(r"domain length  $L(t)$")
    ax.set_title("Coarsening exponents: diffusion-limited $t^{1/3}$ vs\n"
                 "conversion-limited $t^{1/2}$  (Lee 2021)")
    ax.legend(fontsize=9, loc="upper left")
    fig.tight_layout()
    out = FIGURES / fname
    fig.savefig(out)
    plt.close(fig)
    print(f"  saved {out}  (Model B {eB:.3f},  Model A {eA:.3f})")


if __name__ == "__main__":
    import sys
    if "--selftest" in sys.argv:
        for kr in [0.0, 0.02, 0.1, 0.3]:
            sim = FuelDrivenCH2D(N=96, L=96.0, k_react=kr); sim.initialize(seed=3)
            sim.run(2000)
            print(f"k_react={kr}: L={sim.domain_length():.2f}  "
                  f"dense={sim.dense_fraction():.3f}  mean={sim.phi.mean():+.4f}")
        print(f"k_diss (a=kappa=M=1) = {k_dissolution(-1.0,1.0,1.0):.4f}")
        sys.exit()

    print("Fuel-driven droplets : reaction-arrested coarsening")
    print("=" * 52)
    print("\n[1] Rate sweep for the bifurcation diagram (critical quench)")
    sweep = run_rate_sweep([0.01, 0.02, 0.04, 0.07, 0.10, 0.14, 0.18, 0.22,
                            0.25, 0.30], T=1500.0, N=192, L=160.0)
    plot_bifurcation(sweep)

    print("\n[2] Passive vs fuel-driven morphology (droplet quench)")
    sweep_d = run_rate_sweep([0.0, 0.05], phi_mean=-0.35, phi_ss=-0.35,
                             T=1200.0, N=192, L=160.0)
    plot_passive_vs_fuel(sweep_d, L=160.0)

    print("\n[3] Coarsening exponent (Module E)")
    plot_coarsening_exponent()
    print("\nDone.")
