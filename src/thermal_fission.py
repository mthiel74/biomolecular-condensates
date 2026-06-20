#!/usr/bin/env python3
r"""Thermal-gradient droplet fission: membrane-free protocell division
(Ianeselli, Atienza-Sanz, ... Braun & Mast 2022, Nature Chemistry;
"Water cycles and heat fluxes drive ... division of protocell-like droplets").

Phase 2 of the nonlinear-dynamics study.  A temperature gradient across a
condensate makes the demixing tendency position-dependent: chi ~ 1/T, so the
cold side phase-separates more strongly than the warm side.  In the double-well
Cahn-Hilliard model this is a spatially varying well depth a = a(x) (a controls
the depth of the dense minimum, a ~ chi - chi_c):

    mu = -a(x) phi + b phi^3 - kappa lap(phi)
    d phi/dt = M lap(mu)

A droplet straddling the gradient is squeezed where a is small (warm) and fed
where a is large (cold).  With a warm band through the droplet's waist the
interface necks and pinches -- the droplet DIVIDES into two daughters with no
membrane, purely from the thermal field.  This is a candidate prebiotic
division mechanism (heated rock pores on the early Earth).

The local conservation of the Cahn-Hilliard dynamics guarantees the daughters
inherit the parent material (no mass is created), and the gradient -- not
random noise -- selects where the cut is made.

Figures -> figures/.  Run as a script:  python3 src/thermal_fission.py
"""

from __future__ import annotations

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path

from nonlinear_analysis import set_style, ACCENT, ACCENT2, GOLD, INK, DIV_CMAP

FIGURES = Path(__file__).resolve().parent.parent / "figures"
FIGURES.mkdir(exist_ok=True)


class ThermalGradientCH2D:
    r"""Double-well Cahn-Hilliard with a position-dependent well depth a(x,y).

        mu = -a_field * phi + b phi^3 - kappa lap(phi)
        d phi/dt = M lap(mu)

    a_field encodes the thermal map chi(T(x)) ~ 1/T: large a = cold = strong
    demixing; small (or negative) a = warm = mixed.  The constant Eyre
    stabiliser uses A >= max(a_field) for gradient stability.
    """

    def __init__(self, L=128.0, N=256, b=1.0, kappa=1.0, mobility=1.0,
                 dt=0.04, a_field=None, A=None):
        self.L, self.N = L, N
        self.b, self.kappa = b, kappa
        self.M, self.dt = mobility, dt
        self.dx = L / N
        if a_field is None:
            a_field = np.ones((N, N))
        self.a_field = a_field
        self.A = float(np.max(a_field)) + 1.5 if A is None else A
        k = 2 * np.pi * np.fft.fftfreq(N, d=self.dx)
        self.KX, self.KY = np.meshgrid(k, k)
        self.k2 = self.KX**2 + self.KY**2
        self._denom = 1.0 + dt * mobility * self.k2 * (self.A + kappa * self.k2)
        self.phi = None
        self.time = 0.0

    def set_droplet(self, cx, cy, radius, phi_in=1.0, phi_out=-1.0, width=2.0):
        """Initialise a single smooth circular droplet."""
        x = np.linspace(0, self.L, self.N, endpoint=False)
        X, Y = np.meshgrid(x, x)
        r = np.sqrt((X - cx)**2 + (Y - cy)**2)
        self.phi = phi_out + 0.5 * (phi_in - phi_out) * (1 - np.tanh((r - radius) / width))
        self.time = 0.0

    def step(self):
        phi = self.phi
        g = -self.a_field * phi + self.b * phi**3 - self.A * phi
        g_hat = np.fft.fft2(g)
        phi_hat = np.fft.fft2(phi)
        self.phi = np.real(np.fft.ifft2(
            (phi_hat - self.dt * self.M * self.k2 * g_hat) / self._denom))
        self.time += self.dt

    def run(self, n_steps):
        for _ in range(n_steps):
            self.step()

    def droplet_count(self, threshold=0.0):
        from scipy.ndimage import label
        return label(self.phi > threshold)[1]

    def dense_mass(self, threshold=0.0):
        return float(np.sum(self.phi[self.phi > threshold]) * self.dx**2)


# ---------------------------------------------------------------------------
# Thermal maps  a(x,y) ~ chi(T(x)) ~ 1/T
# ---------------------------------------------------------------------------

def warm_band_field(N, L, a_cold=1.0, a_warm=-0.6, y0=None, sigma=10.0):
    """Cold everywhere except a warm horizontal band at y0 (a dips there).

    A droplet centred on the band is fed from above and below (cold) but
    dissolved at its waist (warm) -> it necks and divides."""
    if y0 is None:
        y0 = L / 2
    y = np.linspace(0, L, N, endpoint=False)
    X, Y = np.meshgrid(y, y)
    band = np.exp(-((Y - y0)**2) / (2 * sigma**2))
    a = a_cold - (a_cold - a_warm) * band
    return a


def linear_gradient_field(N, L, a_lo=0.2, a_hi=1.6):
    """Monotonic thermal gradient a increasing with y (cold top, warm bottom)."""
    y = np.linspace(0, L, N, endpoint=False)
    X, Y = np.meshgrid(y, y)
    return a_lo + (a_hi - a_lo) * Y / L


# ---------------------------------------------------------------------------
# Fission run + figures
# ---------------------------------------------------------------------------

# A configuration that divides cleanly (verified): a cold matrix with a warm
# band through the droplet's waist; weak thermal noise breaks the metastable
# neck so the pinch-off completes.
FISSION = dict(a_cold=1.3, a_warm=-3.5, sigma=10.0, radius=28.0,
               dt=0.02, noise=0.03)


def run_fission(N=256, L=128.0, total_steps=16000, record_steps=None,
                cfg=None, seed=0):
    """Evolve a single droplet on a warm-band thermal map; return snapshots."""
    cfg = {**FISSION, **(cfg or {})}
    af = warm_band_field(N, L, a_cold=cfg["a_cold"], a_warm=cfg["a_warm"],
                         sigma=cfg["sigma"])
    sim = ThermalGradientCH2D(N=N, L=L, a_field=af, dt=cfg["dt"])
    sim.set_droplet(L / 2, L / 2, radius=cfg["radius"])
    rng = np.random.default_rng(seed)
    sim.phi = sim.phi + cfg["noise"] * rng.standard_normal((N, N))
    if record_steps is None:
        record_steps = [0, total_steps // 8, total_steps // 3,
                        total_steps // 2, total_steps]
    snaps = []
    prev = 0
    for ns in record_steps:
        sim.run(ns - prev)
        prev = ns
        from scipy.ndimage import label
        snaps.append((sim.time, sim.phi.copy(), label(sim.phi > 0.3)[1]))
    return af, snaps, sim


def plot_fission_sequence(N=256, L=128.0, fname="thermal_fission_sequence.png"):
    """Thermal map + time sequence of a dividing droplet."""
    set_style()
    af, snaps, _ = run_fission(N=N, L=L)
    ncol = len(snaps)
    fig, axes = plt.subplots(1, ncol + 1, figsize=(3.0 * (ncol + 1), 3.4))

    # thermal map a(x) ~ chi ~ 1/T
    im0 = axes[0].imshow(af, cmap="inferno", origin="lower",
                         extent=[0, L, 0, L])
    axes[0].set_title("thermal map\n$a(x)\\sim\\chi\\sim1/T$", fontsize=11)
    axes[0].set_xticks([]); axes[0].set_yticks([])
    cb = fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)
    cb.ax.tick_params(labelsize=7)
    axes[0].text(0.5, 0.5, "warm", color="white", ha="center", fontsize=8,
                 transform=axes[0].transAxes)

    for ax, (t, phi, n) in zip(axes[1:], snaps):
        ax.imshow(phi, cmap=DIV_CMAP, vmin=-1.1, vmax=1.1, origin="lower",
                  extent=[0, L, 0, L])
        ax.set_title(fr"$t={t:.0f}$   ({n} drop{'s' if n != 1 else ''})",
                     fontsize=11)
        ax.set_xticks([]); ax.set_yticks([])

    fig.suptitle("Thermal-gradient fission: a temperature band divides a "
                 "membrane-free droplet  (Ianeselli 2022)",
                 fontsize=14, y=1.04)
    fig.tight_layout()
    out = FIGURES / fname
    fig.savefig(out)
    plt.close(fig)
    print(f"  saved {out}  (final droplet count {snaps[-1][2]})")


def scan_division(N=128, L=96.0, total_steps=11000,
                  warm_depths=None, radii=None, seed=0):
    """For each (warm-band depth, droplet radius) record the daughter count."""
    if warm_depths is None:
        warm_depths = np.array([-0.5, -1.2, -2.0, -2.8, -3.6, -4.4])
    if radii is None:
        radii = np.array([12, 16, 20, 24, 28, 32], dtype=float)
    from scipy.ndimage import label
    grid = np.zeros((len(warm_depths), len(radii)))
    for i, aw in enumerate(warm_depths):
        for j, r in enumerate(radii):
            cfg = dict(a_warm=aw, radius=r)
            _, snaps, sim = run_fission(N=N, L=L, total_steps=total_steps,
                                        record_steps=[total_steps], cfg=cfg,
                                        seed=seed)
            # divided = 2+ daughters each carrying real mass
            n = label(sim.phi > 0.3)[1]
            grid[i, j] = n
        print(f"  warm={aw:+.1f}: counts={grid[i].astype(int).tolist()}",
              flush=True)
    return warm_depths, radii, grid


def plot_division_phase_diagram(fname="thermal_division_phase_diagram.png",
                                cached=None):
    """Phase diagram: gradient strength vs droplet size -> division threshold."""
    set_style()
    if cached is not None:
        warm, radii, grid = cached
    else:
        warm, radii, grid = scan_division()
    fig, ax = plt.subplots(figsize=(7.6, 6.0))
    grad = -warm                      # gradient strength = warm-band depth
    divided = (grid >= 2).astype(float)
    im = ax.imshow(divided, origin="lower", aspect="auto", cmap="RdBu",
                   vmin=0, vmax=1,
                   extent=[radii[0], radii[-1], grad[0], grad[-1]])
    # overlay the actual daughter counts
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            ax.text(radii[j], grad[i], int(grid[i, j]), ha="center",
                    va="center", fontsize=9, color=INK)
    ax.set_xlabel(r"droplet radius  $R$")
    ax.set_ylabel(r"thermal-gradient strength  $|a_{\rm warm}|$")
    ax.set_title("Division threshold: gradient strength vs droplet size\n"
                 "(blue = divides into $\\geq2$ daughters)")
    fig.colorbar(im, ax=ax, ticks=[0, 1], label="divides (1) / intact (0)")
    fig.tight_layout()
    out = FIGURES / fname
    fig.savefig(out)
    plt.close(fig)
    print(f"  saved {out}")


def animate_fission(N=200, L=128.0, total_steps=16000, n_frames=80,
                    fname="thermal_fission.mp4"):
    """Animation of the division process."""
    import matplotlib.animation as animation
    set_style()
    cfg = FISSION
    af = warm_band_field(N, L, a_cold=cfg["a_cold"], a_warm=cfg["a_warm"],
                         sigma=cfg["sigma"])
    sim = ThermalGradientCH2D(N=N, L=L, a_field=af, dt=cfg["dt"])
    sim.set_droplet(L / 2, L / 2, radius=cfg["radius"])
    rng = np.random.default_rng(0)
    sim.phi = sim.phi + cfg["noise"] * rng.standard_normal((N, N))

    fig, ax = plt.subplots(figsize=(5.2, 5.2))
    im = ax.imshow(sim.phi, cmap=DIV_CMAP, vmin=-1.1, vmax=1.1, origin="lower",
                   extent=[0, L, 0, L])
    ax.set_xticks([]); ax.set_yticks([])
    ttl = ax.set_title("")
    steps_per = max(1, total_steps // n_frames)
    from scipy.ndimage import label

    def update(f):
        sim.run(steps_per)
        im.set_data(sim.phi)
        ttl.set_text(f"thermal fission   $t={sim.time:.0f}$   "
                     f"({label(sim.phi > 0.3)[1]} droplets)")
        return im, ttl

    anim = animation.FuncAnimation(fig, update, frames=n_frames, blit=False)
    out = FIGURES / fname
    try:
        anim.save(out, writer="ffmpeg", fps=18, dpi=110)
        print(f"  saved {out}")
    except Exception as e:
        print(f"  [skip animation: {e}]")
    plt.close(fig)


if __name__ == "__main__":
    print("Thermal-gradient droplet fission")
    print("=" * 40)
    print("\n[1] Division sequence")
    plot_fission_sequence()
    print("\n[2] Division phase diagram")
    plot_division_phase_diagram()
    print("\n[3] Animation")
    animate_fission()
    print("\nDone.")
