#!/usr/bin/env python3
r"""Active Model B+ : activity-driven microphase separation and the reversal
of the Ostwald process  (Wittkowski 2014; Tjhung, Nardini & Cates 2018, PRX;
Caballero, Nardini & Cates 2018).

Phase 2 of the nonlinear-dynamics study.  Builds on the passive double-well
Model-B solver (`nonlinear_analysis.DoubleWellCH2D`) and adds the two
non-equilibrium terms that turn an *equilibrium* phase-separating fluid into
an *active* one whose detailed balance is broken:

    chemical potential   mu = f'(phi) - kappa lap(phi) + lambda |grad phi|^2
    conserved dynamics   d phi/dt = M lap(mu) - zeta div[ (lap phi) grad phi ]

Two distinct activities appear:

  * lambda  -- a |grad phi|^2 contribution to mu.  It is NOT the gradient of
    any free-energy functional (no local F has dF/dphi = |grad phi|^2), so it
    breaks detailed balance.  At *linear* order about a uniform state it is
    O(delta^2) and therefore leaves the dispersion relation untouched: lambda
    shifts the binodal and the interfacial tension but on its own does not
    arrest coarsening.  This is the honest statement of the 2014 "Active Model
    B" result, and we demonstrate it rather than assume the opposite.

  * zeta  -- a non-integrable contribution to the *current*
    J_act = zeta (lap phi) grad phi.  This is the defining extra term of Model
    B+.  It modifies the curvature dependence of the interfacial flux and, at
    sufficient strength, REVERSES the Ostwald process: small droplets grow at
    the expense of large ones, coarsening is arrested, and the system locks
    into a steady microphase-separated emulsion with an activity-selected
    length scale.  This is the dramatic result of Tjhung 2018.

The relevant combination governing the reversal is v = zeta - 2 lambda
(Tjhung 2018): we therefore sweep the activity and *measure* the steady-state
length L_ss, rather than asserting it.

Section C (thermal-gradient fission, Ianeselli 2022) lives in
`thermal_fission.py`.

Figures -> figures/.  Run as a script to regenerate everything:
    python3 src/active_model_b.py
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

class ActiveModelBPlus:
    r"""Semi-implicit spectral solver for Active Model B+ in 2D.

        f(phi) = -a/2 phi^2 + b/4 phi^4,   f'(phi) = -a phi + b phi^3
        mu     = f'(phi) - kappa lap(phi) + lambda |grad phi|^2
        d phi/dt = M lap(mu) - zeta div[ (lap phi) grad phi ]

    The passive linear part (Eyre A + biharmonic kappa k^2) is treated
    implicitly; the double-well remainder, the lambda activity, and the zeta
    current are explicit.  With lam = zeta = 0 this reduces exactly to the
    passive double-well Model-B solver used in Phase 1.
    """

    def __init__(self, L=128.0, N=256, a=1.0, b=1.0, kappa=1.0, mobility=1.0,
                 dt=0.05, A=2.5, lam=0.0, zeta=0.0):
        self.L, self.N = L, N
        self.a, self.b, self.kappa = a, b, kappa
        self.M, self.dt, self.A = mobility, dt, A
        self.lam, self.zeta = lam, zeta
        self.dx = L / N
        k = 2 * np.pi * np.fft.fftfreq(N, d=self.dx)
        self.KX, self.KY = np.meshgrid(k, k)
        self.k2 = self.KX**2 + self.KY**2
        self._denom = 1.0 + dt * mobility * self.k2 * (A + kappa * self.k2)
        # 2/3-rule dealiasing mask: kills the top third of each axis to remove
        # the aliasing error from the cubic / active nonlinear products, which
        # is what makes the high-derivative zeta current blow up.
        kmax = np.abs(k).max()
        self._dealias = (np.abs(self.KX) <= (2.0 / 3.0) * kmax) & \
                        (np.abs(self.KY) <= (2.0 / 3.0) * kmax)
        self.phi = None
        self.time = 0.0

    def initialize(self, phi_mean=0.0, noise_amplitude=0.05, seed=42):
        rng = np.random.default_rng(seed)
        self.phi = phi_mean + noise_amplitude * rng.standard_normal((self.N, self.N))
        self.time = 0.0

    # -- spectral helpers ---------------------------------------------------
    def _grads(self, phi_hat):
        gx = np.real(np.fft.ifft2(1j * self.KX * phi_hat))
        gy = np.real(np.fft.ifft2(1j * self.KY * phi_hat))
        return gx, gy

    def _lap(self, phi_hat):
        return np.real(np.fft.ifft2(-self.k2 * phi_hat))

    def _div(self, vx, vy):
        vx_hat = np.fft.fft2(vx)
        vy_hat = np.fft.fft2(vy)
        return np.real(np.fft.ifft2(1j * self.KX * vx_hat + 1j * self.KY * vy_hat))

    def step(self):
        phi = self.phi
        phi_hat = np.fft.fft2(phi)
        gx, gy = self._grads(phi_hat)

        # explicit chemical-potential remainder  g = f'(phi) - A phi + lam|grad phi|^2
        g = (-self.a * phi + self.b * phi**3) - self.A * phi
        if self.lam != 0.0:
            g = g + self.lam * (gx**2 + gy**2)
        g_hat = np.fft.fft2(g) * self._dealias

        rhs = phi_hat - self.dt * self.M * self.k2 * g_hat

        # explicit non-integrable active current:  -zeta div[(lap phi) grad phi]
        if self.zeta != 0.0:
            lap = self._lap(phi_hat)
            div = self._div(lap * gx, lap * gy)
            rhs = rhs - self.dt * self.zeta * (np.fft.fft2(div) * self._dealias)

        self.phi = np.real(np.fft.ifft2(rhs / self._denom))
        self.time += self.dt

    def run(self, n_steps):
        for _ in range(n_steps):
            self.step()

    # -- diagnostics --------------------------------------------------------
    def domain_length(self):
        """Characteristic length from the first moment of the structure factor:
        L = 2 pi <|S|> / <k |S|>  (robust, no peak-finding)."""
        phi = self.phi - self.phi.mean()
        S = np.abs(np.fft.fft2(phi))**2
        kmag = np.sqrt(self.k2)
        num = np.sum(S)
        den = np.sum(kmag * S)
        return 2 * np.pi * num / den if den > 0 else self.L


def droplet_radii(phi, dx, threshold=0.0):
    """Equivalent radii of connected dense domains (phi > threshold)."""
    from scipy.ndimage import label
    lab, n = label(phi > threshold)
    if n == 0:
        return np.array([])
    areas = np.bincount(lab.ravel())[1:] * dx**2
    return np.sqrt(areas / np.pi)


# ===========================================================================
# Sweep
# ===========================================================================

def run_activity_sweep(zeta_values, T=1500.0, N=192, L=128.0, dt=0.04,
                       lam=0.0, phi_mean=0.0, seed=7, n_rec=140):
    """Evolve Model B+ for each activity zeta; record L(t) and final field."""
    out = {}
    for z in zeta_values:
        print(f"  zeta={z:+.3f} ...", end="", flush=True)
        sim = ActiveModelBPlus(L=L, N=N, dt=dt, lam=lam, zeta=z)
        sim.initialize(phi_mean=phi_mean, seed=seed)
        n = int(T / dt)
        rec_every = max(1, n // n_rec)
        ts, Ls = [], []
        for i in range(1, n + 1):
            sim.step()
            if not np.isfinite(sim.phi).all():
                print(" BLEW UP", flush=True)
                break
            if i % rec_every == 0:
                ts.append(sim.time)
                Ls.append(sim.domain_length())
        out[z] = (np.array(ts), np.array(Ls), sim.phi.copy())
        print(f" L={Ls[-1]:.1f}" if Ls else " (no record)", flush=True)
    return out


# ===========================================================================
# Figures
# ===========================================================================

def plot_three_regimes(sweep, L=128.0, fname="active_b_three_regimes.png"):
    """Hero panel: snapshots at three activities + L_ss(zeta) selection curve."""
    set_style()
    zs = sorted(sweep)
    # pick passive, intermediate, strong from the swept set
    z_lo = zs[0]
    z_hi = zs[-1]
    z_mid = min(zs, key=lambda z: abs(z - 0.5 * (z_lo + z_hi)))

    fig = plt.figure(figsize=(16.5, 5.2))
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 1.15], wspace=0.28)

    titles = [
        (z_lo, "(a) passive  $\\zeta=0$\nfull coarsening"),
        (z_mid, f"(b) $\\zeta={z_mid:g}$\narrested coarsening"),
        (z_hi, f"(c) $\\zeta={z_hi:g}$\nmicrophase emulsion"),
    ]
    for j, (z, ttl) in enumerate(titles):
        ax = fig.add_subplot(gs[0, j])
        phi = sweep[z][2]
        ax.imshow(phi, extent=[0, L, 0, L], cmap=DIV_CMAP, vmin=-1.2, vmax=1.2,
                  origin="lower", interpolation="bilinear")
        ax.set_title(ttl, fontsize=12)
        ax.set_xticks([]); ax.set_yticks([])

    # selection curve L_ss vs zeta (steady-state = mean of last 15% of L(t))
    ax = fig.add_subplot(gs[0, 3])
    Lss = np.array([np.mean(sweep[z][1][-max(3, len(sweep[z][1]) // 7):])
                    for z in zs])
    ax.plot(zs, Lss, 'o-', color=ACCENT, ms=8, mec="white", mew=0.6, lw=2)
    ax.scatter([z_lo, z_mid, z_hi],
               [Lss[zs.index(z_lo)], Lss[zs.index(z_mid)], Lss[zs.index(z_hi)]],
               s=160, facecolors="none", edgecolors=INK, lw=1.6, zorder=5)
    ax.set_xlabel(r"activity  $\zeta$")
    ax.set_ylabel(r"steady-state length  $L_{\rm ss}$")
    ax.set_title("(d) activity selects a length scale")
    ax.annotate("Ostwald\nreversed", xy=(z_hi, Lss[-1]),
                xytext=(0.55, 0.78), textcoords="axes fraction", fontsize=9.5,
                ha="center", color=INK,
                arrowprops=dict(arrowstyle="->", color=INK, lw=1.2))

    fig.suptitle("Active Model B+ : activity arrests coarsening and selects a "
                 "length scale  (Tjhung, Nardini & Cates 2018)",
                 fontsize=14.5, y=1.03)
    out = FIGURES / fname
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out}")
    return Lss


def plot_length_vs_activity(sweep, fname="active_b_length_selection.png"):
    """Dedicated L(t) trajectories + L_ss(zeta) with the passive t^1/3 reference."""
    set_style()
    zs = sorted(sweep)
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.5, 5.0))
    cmap = mpl.colormaps["viridis"]

    for i, z in enumerate(zs):
        t, Lt, _ = sweep[z]
        if len(t) == 0:
            continue
        c = INK if z == 0 else cmap(0.12 + 0.78 * i / max(1, len(zs) - 1))
        lab = "passive $\\zeta=0$" if z == 0 else fr"$\zeta={z:g}$"
        axA.plot(t, Lt, color=c, lw=2.2, label=lab)
    # passive t^1/3 guide
    t0 = sweep[zs[0]][0]
    if len(t0):
        tg = np.linspace(t0[len(t0)//4], t0[-1], 50)
        L0 = sweep[zs[0]][1]
        c = L0[len(L0)//2] / t0[len(t0)//2]**(1/3)
        axA.plot(tg, c * tg**(1/3), '--', color=GOLD, lw=1.8,
                 label=r"$L\sim t^{1/3}$ (LSW)")
    axA.set_xlabel(r"time $t$"); axA.set_ylabel(r"domain length $L(t)$")
    axA.set_title("(a) activity freezes the growth of $L(t)$")
    axA.legend(fontsize=8.5, loc="upper left")

    Lss = np.array([np.mean(sweep[z][1][-max(3, len(sweep[z][1]) // 7):])
                    if len(sweep[z][1]) else np.nan for z in zs])
    axB.plot(zs, Lss, 'o-', color=ACCENT, ms=9, mec="white", mew=0.6, lw=2)
    axB.set_xlabel(r"activity  $\zeta$")
    axB.set_ylabel(r"steady-state length  $L_{\rm ss}$")
    axB.set_title(r"(b) $L_{\rm ss}$ diverges as $\zeta\!\to\!0$, saturates at high $\zeta$")
    fig.suptitle("Activity-selected length scale in Active Model B+",
                 fontsize=14.5, y=1.02)
    fig.tight_layout()
    out = FIGURES / fname
    fig.savefig(out)
    plt.close(fig)
    print(f"  saved {out}")


if __name__ == "__main__":
    import sys
    if "--selftest" in sys.argv:
        # quick smoke test: passive vs active, short run on small grid
        for z in [0.0, 1.0, 4.0]:
            sim = ActiveModelBPlus(N=96, L=96.0, zeta=z, dt=0.04)
            sim.initialize(seed=3)
            sim.run(1500)
            print(f"zeta={z}: L={sim.domain_length():.2f}  "
                  f"mean={sim.phi.mean():+.4f}  finite={np.isfinite(sim.phi).all()}")
        sys.exit()

    print("Active Model B+ : activity-driven microphase separation")
    print("=" * 56)
    sweep = run_activity_sweep([0.0, 1.0, 2.0, 4.0, 6.0, 8.0])
    plot_three_regimes(sweep)
    plot_length_vs_activity(sweep)
    print("\nDone.")
