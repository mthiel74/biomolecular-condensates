#!/usr/bin/env python3
r"""Wet-dry cycling: prebiotic tidal pools stabilise a population of droplets
(Ianeselli 2022, Nat. Chem.; Campbell/Damer & Deamer 2020 wet-dry hypothesis;
Ross & Deamer 2016).

Phase 3 of the nonlinear-dynamics study.  On the early Earth, rock pools and
mineral surfaces went through repeated WET (flooded, dilute) and DRY
(evaporated, concentrated) cycles.  In the Flory-Huggins picture the effective
interaction tracks solvent activity / temperature, so cycling drives

    chi(t) = chi_0 + d_chi * sin(2 pi t / T) .

In the double-well Cahn-Hilliard model chi maps onto the well depth
a ~ chi - chi_c, so a(t) swings the system between:

  * WET  (a < 0): one-phase well -- condensates DISSOLVE, contents remix;
  * DRY  (a > 0): double well   -- the (off-critical) mixture is spinodally
                  unstable and RE-NUCLEATES a fresh crop of many droplets.

The consequence is striking.  Held at constant (dry) chi the droplets undergo
Ostwald ripening: big droplets eat small ones and the population collapses
toward a single winner (droplet number ~ t^{-2/3}).  Under wet-dry cycling each
dissolution ERASES the size differences that drive ripening; the next dry phase
re-nucleates many droplets before coarsening can run away.  A POPULATION of
droplets is maintained indefinitely -- and a passive cargo is re-concentrated
on every dry phase (a molecular ratchet for prebiotic chemistry).

Figures -> figures/.  Run as a script:
    python3 src/wet_dry_cycling.py
    python3 src/wet_dry_cycling.py --quick
"""

from __future__ import annotations

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.ndimage import label

from nonlinear_analysis import set_style, ACCENT, ACCENT2, GOLD, INK, DIV_CMAP

FIGURES = Path(__file__).resolve().parent.parent / "figures"
FIGURES.mkdir(exist_ok=True)


# ===========================================================================
# Cahn-Hilliard solver with time-dependent well depth a(t) + passive cargo
# ===========================================================================

class CyclingCH2D:
    r"""Double-well Model-B with a(t), an optional passive cargo c that
    partitions into the dense phase (drift up grad phi).

        mu = -a(t) phi + b phi^3 - kappa lap(phi)
        d phi/dt = M lap(mu)
        d c/dt   = D_c lap(c) - alpha D_c div(c grad phi)
    """

    def __init__(self, L=160.0, N=128, b=1.0, kappa=1.0, mobility=1.0,
                 dt=0.06, A=4.6, D_c=0.5, alpha=2.4, noise=0.0015, seed=0):
        self.L, self.N, self.b, self.kappa = L, N, b, kappa
        self.M, self.dt, self.A = mobility, dt, A
        self.D_c, self.alpha, self.noise = D_c, alpha, noise
        self.dx = L / N
        k = 2 * np.pi * np.fft.fftfreq(N, d=self.dx)
        self.KX, self.KY = np.meshgrid(k, k)
        self.k2 = self.KX**2 + self.KY**2
        self._denom = 1.0 + dt * mobility * self.k2 * (A + kappa * self.k2)
        self._diff_if = np.exp(-D_c * self.k2 * dt)
        self.rng = np.random.default_rng(seed)
        self.phi = None
        self.c = None
        self.time = 0.0

    def initialize(self, phi_mean=-0.40, c_mean=0.30, amp=0.05):
        self.phi = phi_mean + amp * self.rng.standard_normal((self.N, self.N))
        self.c = np.full((self.N, self.N), c_mean)
        self.time = 0.0

    def _grad(self, u):
        uh = np.fft.fft2(u)
        return (np.real(np.fft.ifft2(1j * self.KX * uh)),
                np.real(np.fft.ifft2(1j * self.KY * uh)))

    def _lap(self, u):
        return np.real(np.fft.ifft2(-self.k2 * np.fft.fft2(u)))

    def step(self, a):
        phi = self.phi
        g = (-a * phi + self.b * phi**3) - self.A * phi
        phi_hat = np.fft.fft2(phi)
        phi = np.real(np.fft.ifft2(
            (phi_hat - self.dt * self.M * self.k2 * np.fft.fft2(g))
            / self._denom))
        if self.noise:
            phi = phi + self.noise * self.rng.standard_normal((self.N, self.N))
        # cargo transport (drift into dense phase)
        if self.c is not None:
            gx, gy = self._grad(self.phi)
            cx, cy = self._grad(self.c)
            drift = -self.alpha * self.D_c * (cx * gx + cy * gy
                                              + self.c * self._lap(self.phi))
            c_hat = self._diff_if * (np.fft.fft2(self.c)
                                     + self.dt * np.fft.fft2(drift))
            self.c = np.maximum(np.real(np.fft.ifft2(c_hat)), 0.0)
        self.phi = phi
        self.time += self.dt

    def droplet_count(self, thresh=0.0, min_px=6):
        lab, n = label(self.phi > thresh)
        if n == 0:
            return 0
        sizes = np.bincount(lab.ravel())[1:]
        return int(np.sum(sizes >= min_px))

    def cargo_enrichment(self):
        """Peak cargo concentration relative to the mean (fold enrichment)."""
        return float(self.c.max() / max(self.c.mean(), 1e-9))


def chi_wave(t, a0, d_a, period):
    return a0 + d_a * np.sin(2 * np.pi * t / period)


# ===========================================================================
# Runs
# ===========================================================================

def run_protocol(cycling, a0=1.0, d_a=1.8, period=300.0, T_total=4200.0,
                 N=128, L=160.0, dt=0.06, n_rec=240, seed=0, snap_times=None):
    """Evolve either constant (a=a0) or cycling a(t); record droplet count,
    cargo enrichment, a(t), and a handful of field snapshots."""
    sim = CyclingCH2D(N=N, L=L, dt=dt, seed=seed)
    sim.initialize()
    n = int(T_total / dt)
    rec_every = max(1, n // n_rec)
    if snap_times is None:
        snap_times = np.linspace(0.15, 1.0, 5) * T_total
    ts, counts, enrich, avals = [], [], [], []
    snaps, snap_ts = [], []
    si = 0
    for i in range(1, n + 1):
        a = a0 if not cycling else chi_wave(sim.time, a0, d_a, period)
        sim.step(a)
        if i % rec_every == 0:
            ts.append(sim.time)
            counts.append(sim.droplet_count())
            enrich.append(sim.cargo_enrichment())
            avals.append(a)
        if si < len(snap_times) and sim.time >= snap_times[si]:
            snaps.append(sim.phi.copy()); snap_ts.append(sim.time); si += 1
    return dict(t=np.array(ts), count=np.array(counts),
                enrich=np.array(enrich), a=np.array(avals),
                snaps=snaps, snap_ts=snap_ts, cycling=cycling,
                a0=a0, d_a=d_a, period=period, L=L)


# ===========================================================================
# Figures
# ===========================================================================

def plot_population_metric(const, cyc, fname="wetdry_population.png"):
    """Hero: droplet number vs time (constant collapses to one; cycling stays
    a population), the chi(t) forcing, and cargo re-concentration per cycle."""
    set_style()
    fig = plt.figure(figsize=(13.5, 7.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.85], hspace=0.34,
                          wspace=0.24)
    axN = fig.add_subplot(gs[0, :])
    axW = fig.add_subplot(gs[1, 0])
    axC = fig.add_subplot(gs[1, 1])

    # (a) droplet number
    axN.plot(const["t"], const["count"], '-', color=ACCENT2, lw=2.4,
             label="constant $\\chi$ (Ostwald ripening)")
    axN.plot(cyc["t"], cyc["count"], '-', color=ACCENT, lw=2.0,
             label="wet--dry cycling")
    # t^{-2/3} guide for the constant (coarsening) curve
    m = const["count"] > 0
    tg = const["t"][m]
    if len(tg) > 10:
        c0 = const["count"][m][len(tg) // 6] * (tg[len(tg) // 6])**(2 / 3)
        axN.plot(tg, c0 * tg**(-2 / 3), ':', color=INK, lw=1.6,
                 label=r"$N_{\rm drop}\sim t^{-2/3}$ (LSW)")
    # shade dry/wet half-cycles
    for k in range(int(cyc["t"][-1] / cyc["period"]) + 1):
        t1 = k * cyc["period"]
        axN.axvspan(t1, t1 + cyc["period"] / 2, color=GOLD, alpha=0.06)
    axN.set_xlabel("time  $t$")
    axN.set_ylabel("number of droplets")
    axN.set_title("(a) wet--dry cycling sustains a droplet population; "
                  "constant $\\chi$ coarsens toward one winner")
    axN.legend(loc="upper right", fontsize=10)
    axN.set_ylim(0, max(cyc["count"].max(), 5) * 1.12)

    # (b) the forcing waveform
    axW.plot(cyc["t"], cyc["a"], '-', color=ACCENT, lw=2.0)
    axW.axhline(0, color=INK, lw=1.0, ls="--")
    axW.fill_between(cyc["t"], cyc["a"], 0, where=cyc["a"] > 0, color=ACCENT,
                     alpha=0.12)
    axW.fill_between(cyc["t"], cyc["a"], 0, where=cyc["a"] < 0, color=ACCENT2,
                     alpha=0.12)
    axW.text(cyc["period"] * 0.25, cyc["d_a"] * 0.6, "DRY\n(condense)",
             ha="center", fontsize=9, color="#9a2a4a")
    axW.text(cyc["period"] * 0.75, -cyc["d_a"] * 0.6, "WET\n(dissolve)",
             ha="center", fontsize=9, color=ACCENT2)
    axW.set_xlim(0, 2.2 * cyc["period"])
    axW.set_xlabel("time  $t$")
    axW.set_ylabel(r"well depth  $a\sim\chi-\chi_c$")
    axW.set_title(r"(b) forcing  $\chi(t)=\chi_0+\Delta\chi\,\sin(2\pi t/T)$")

    # (c) cargo enrichment ratchet
    axC.plot(const["t"], const["enrich"], '-', color=ACCENT2, lw=2.0,
             label="constant $\\chi$")
    axC.plot(cyc["t"], cyc["enrich"], '-', color=ACCENT, lw=2.0,
             label="cycling (re-concentrates)")
    axC.set_xlabel("time  $t$")
    axC.set_ylabel(r"cargo enrichment  $c_{\max}/\langle c\rangle$")
    axC.set_title("(c) each dry phase re-concentrates cargo")
    axC.legend(fontsize=9, loc="upper left")

    fig.suptitle("Wet--dry cycling stabilises a population of protocell "
                 "droplets  (Ianeselli 2022; prebiotic tidal pools)",
                 fontsize=14.5, y=0.995)
    out = FIGURES / fname
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out}")


def plot_morphology_comparison(const, cyc, fname="wetdry_morphology.png"):
    """Snapshot rows: constant chi coarsening to one droplet vs cycling keeping
    many."""
    set_style()
    ncol = min(len(const["snaps"]), len(cyc["snaps"]))
    fig, axes = plt.subplots(2, ncol, figsize=(3.0 * ncol, 6.3))
    L = const["L"]
    for j in range(ncol):
        for row, run, lab in [(0, const, "constant $\\chi$"),
                              (1, cyc, "wet--dry cycling")]:
            ax = axes[row, j]
            ax.imshow(run["snaps"][j], cmap=DIV_CMAP, vmin=-1.1, vmax=1.1,
                      origin="lower", extent=[0, L, 0, L],
                      interpolation="bilinear")
            ax.set_xticks([]); ax.set_yticks([])
            if row == 0:
                ax.set_title(fr"$t={run['snap_ts'][j]:.0f}$", fontsize=11)
            if j == 0:
                ax.set_ylabel(lab, fontsize=12)
    axes[0, 0].annotate("coarsens to one winner", (0.5, 1.32),
                        xycoords="axes fraction", ha="center", fontsize=11,
                        color=ACCENT2)
    axes[1, 0].annotate("many droplets persist", (0.5, 1.06),
                        xycoords="axes fraction", ha="center", fontsize=11,
                        color=ACCENT)
    fig.suptitle("Constant vs cycling interaction: Ostwald ripening "
                 "(top) is defeated by periodic dissolution (bottom)",
                 fontsize=13.5, y=1.0)
    fig.tight_layout()
    out = FIGURES / fname
    fig.savefig(out)
    plt.close(fig)
    print(f"  saved {out}")


if __name__ == "__main__":
    import sys
    quick = "--quick" in sys.argv
    kw = dict(N=104, L=128.0, T_total=2400.0, period=300.0, n_rec=180) if quick \
        else dict(N=140, L=150.0, T_total=4500.0, period=360.0, n_rec=260)

    # snapshots at DRY peaks t=(k+0.25)T so the cycling row always shows droplets
    period = kw["period"]
    peaks = np.array([(k + 0.25) * period
                      for k in range(int(kw["T_total"] / period))])
    snap_t = peaks[np.linspace(0.6, len(peaks) - 1, 5).astype(int)]

    print("Wet-dry cycling")
    print("=" * 40)
    print("\n[1] Constant-chi run (Ostwald ripening)")
    const = run_protocol(cycling=False, seed=1, snap_times=snap_t, **kw)
    print(f"    droplets: peak {const['count'].max()} -> {const['count'][-1]}")
    print("[2] Cycling-chi run (population maintained)")
    cyc = run_protocol(cycling=True, seed=1, snap_times=snap_t, **kw)
    print(f"    droplets: {cyc['count'][2]} -> {cyc['count'][-1]}  "
          f"(mean last third {cyc['count'][len(cyc['count'])*2//3:].mean():.1f})")

    print("\n[3] Population metric figure")
    plot_population_metric(const, cyc)
    print("[4] Morphology comparison")
    plot_morphology_comparison(const, cyc)
    print("\nDone.")
