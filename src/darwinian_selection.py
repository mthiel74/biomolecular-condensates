#!/usr/bin/env python3
r"""Darwinian selection among protocells -- natural selection from pure physics.

Phase 3 of the nonlinear-dynamics study; the intellectual climax of the
project.  Three coupled physical fields, plus a shared nutrient, reproduce the
Darwinian cycle of variation, heredity, competition and selection WITHOUT any
biology hand-coded in -- selection is an emergent consequence of the dynamics.

The actors
----------
phi(x,t)   condensate phase field (dilute -1, dense +1).  Fuel-driven
           Cahn-Hilliard so that a *population* of droplets is maintained
           against Ostwald ripening (Zwicker 2015), with a coupling to the
           cargo that lets cargo-rich droplets grow:

               mu = -a phi + b phi^3 - kappa lap(phi) - beta c
               d phi/dt = M lap(mu) - k_react (phi - phi_ss)

           The -beta c term lowers the free energy of the dense phase wherever
           replicators accumulate: a droplet that concentrates cargo draws in
           material and GROWS (osmotic coupling).

c(x,t)     replicator concentration.  Diffuses, partitions preferentially into
           the dense phase, and replicates:

               d c/dt = div[ D_c ( grad c - alpha c grad phi ) ]
                        + r0 * g * c * (1 - c/K) * Nutr

           The drift term -alpha c grad phi is a Smoluchowski flux up the phi
           gradient: at equilibrium c ~ exp(alpha phi), i.e. enriched inside
           droplets by exp(2 alpha) -- the concentration that makes prebiotic
           chemistry fast (earlier modules).

g(x,t)     the GENOTYPE: the intrinsic replication efficiency of the replicator
           at x (the heritable trait).  Carried faithfully by the cargo
           (tracked through h = c*g, transported identically to c), with weak
           mutation (genotype diffusion + noise).  Replication makes copies of
           the same g -- heredity.

Nutr(t)    a shared, global nutrient/fuel pool (chemostat).  Replication
           consumes it; it is slowly resupplied.  Because the resource is
           SHARED, the genotype with the largest g out-replicates the rest and
           drives them extinct -- Malthusian competition / competitive
           exclusion (Eigen 1971).  This is the selective force.

The cycle
---------
Every "generation" the world is evolved, then hit with a thermal pulse (warm
stripes, chi ~ 1/T; Ianeselli 2022) that necks and DIVIDES the large droplets
-- protocell fission.  Lineages are tracked by centroid overlap: each daughter
is assigned its parent, droplets that lose their material dissolve (death),
droplets that pinch in two divide.  The cargo-weighted mean genotype -- the
population mean fitness -- is recorded each generation.

Result: mean fitness rises over generations and saturates at a
mutation-selection balance.  Evolution by natural selection, out of three
PDEs and a conservation law.

Figures + animation -> figures/.  Run as a script:
    python3 src/darwinian_selection.py
    python3 src/darwinian_selection.py --quick     # fast, low-res smoke test
"""

from __future__ import annotations

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from pathlib import Path
from scipy.ndimage import label, center_of_mass

from nonlinear_analysis import set_style, ACCENT, ACCENT2, GOLD, INK, DIV_CMAP

FIGURES = Path(__file__).resolve().parent.parent / "figures"
FIGURES.mkdir(exist_ok=True)

# genotype (replication-efficiency) bounds
G_MIN, G_MAX = 0.45, 1.85
FIT_CMAP = "viridis"          # low fitness -> high fitness


# ===========================================================================
# The coupled world
# ===========================================================================

class ProtocellWorld:
    r"""Coupled phase / cargo / genotype solver with a shared nutrient pool.

    phi : spectral semi-implicit Cahn-Hilliard (Eyre split) + fuel reaction
          + cargo osmotic coupling.
    c,h : integrating-factor spectral diffusion with explicit drift + reaction
          (h = c*g carries the genotype).
    Nutr: global scalar, ODE coupled to the total replication flux.
    """

    def __init__(self, L=128.0, N=128, a=1.0, b=1.0, kappa=1.0, mobility=1.0,
                 dt=0.04, A=2.6, k_react=0.02, phi_ss=-0.42, beta=0.30,
                 D_c=0.6, alpha=2.6, r0=0.60, K=5.0, mut=0.012, death=0.16,
                 nutr_supply=0.22, nutr_use=0.05, seed=0):
        self.L, self.N = L, N
        self.a, self.b, self.kappa = a, b, kappa
        self.M, self.dt, self.A = mobility, dt, A
        self.k_react, self.phi_ss, self.beta = k_react, phi_ss, beta
        self.D_c, self.alpha, self.r0, self.K = D_c, alpha, r0, K
        self.mut, self.death = mut, death
        self.nutr_supply, self.nutr_use = nutr_supply, nutr_use
        self.dx = L / N
        k = 2 * np.pi * np.fft.fftfreq(N, d=self.dx)
        self.KX, self.KY = np.meshgrid(k, k)
        self.k2 = self.KX**2 + self.KY**2
        self._denom = 1.0 + dt * (mobility * self.k2 * (A + kappa * self.k2)
                                  + k_react)
        self._diff_if = np.exp(-D_c * self.k2 * dt)      # integrating factor
        self.rng = np.random.default_rng(seed)
        self.phi = None
        self.c = None
        self.h = None           # h = c * g
        self.Nutr = 1.0
        self.time = 0.0

    # -- spectral helpers ---------------------------------------------------
    def _grad(self, u):
        uh = np.fft.fft2(u)
        return (np.real(np.fft.ifft2(1j * self.KX * uh)),
                np.real(np.fft.ifft2(1j * self.KY * uh)))

    def _lap(self, u):
        return np.real(np.fft.ifft2(-self.k2 * np.fft.fft2(u)))

    # -- initial conditions -------------------------------------------------
    def seed_population(self, n_drops=9, radius=10.0, c0=0.6, margin=16.0):
        """Nucleate n_drops droplets at random non-overlapping centres, each
        loaded with cargo of a random genotype (the initial variation)."""
        N, L = self.N, self.L
        x = np.linspace(0, L, N, endpoint=False)
        X, Y = np.meshgrid(x, x)
        self.phi = np.full((N, N), -1.0)
        self.c = np.full((N, N), 1e-3)
        g0 = np.full((N, N), 0.5 * (G_MIN + G_MAX))
        centres = []
        tries = 0
        while len(centres) < n_drops and tries < 4000:
            tries += 1
            cx, cy = self.rng.uniform(margin, L - margin, 2)
            if all((cx - px)**2 + (cy - py)**2 > (2.6 * radius)**2
                   for px, py, _ in centres):
                g = self.rng.uniform(G_MIN + 0.1, G_MAX - 0.1)
                centres.append((cx, cy, g))
        for (cx, cy, g) in centres:
            r = np.sqrt((X - cx)**2 + (Y - cy)**2)
            blob = 0.5 * (1 - np.tanh((r - radius) / 2.0))
            self.phi = np.maximum(self.phi, -1.0 + 2.0 * blob)
            self.c = self.c + c0 * blob
            g0 = np.where(blob > 0.4, g, g0)
        self.phi += 0.02 * self.rng.standard_normal((N, N))
        self.h = self.c * g0
        self.Nutr = 1.0
        self.time = 0.0
        return centres

    # -- one explicit/semi-implicit time step -------------------------------
    def step(self, replicate=True, a_field=None):
        phi, c, h = self.phi, self.c, self.h
        a = self.a if a_field is None else a_field

        # ---- phase field: semi-implicit CH + fuel + cargo coupling --------
        g_dw = (-a * phi + self.b * phi**3) - self.beta * c - self.A * phi
        phi_hat = np.fft.fft2(phi)
        rhs = phi_hat - self.dt * self.M * self.k2 * np.fft.fft2(g_dw)
        if self.k_react != 0.0:
            rhs[0, 0] += self.dt * self.k_react * self.phi_ss * self.N**2
        phi_new = np.real(np.fft.ifft2(rhs / self._denom))

        # ---- cargo + genotype transport (drift up grad phi) ---------------
        gx, gy = self._grad(phi)
        lap_phi = self._lap(phi)

        def transport(u):
            ux, uy = self._grad(u)
            # -div[ alpha u grad phi ] = -alpha ( grad u . grad phi + u lap phi )
            drift = -self.alpha * (ux * gx + uy * gy + u * lap_phi)
            return self.D_c * drift            # diffusion handled by IF below

        gross = np.zeros_like(c)
        g = self._genotype()
        if replicate:
            growth = np.clip(1.0 - c / self.K, 0.0, 1.0)
            gross = self.r0 * g * c * growth * self.Nutr      # fuelled replication
        # chemostat dilution removes cargo uniformly: genotypes whose growth
        # rate r0*g*Nutr falls below the death rate are out-competed (extinct).
        Ec = transport(c) + gross - self.death * c
        Eh = transport(h) + g * gross - self.death * h        # offspring carry g

        c_hat = self._diff_if * (np.fft.fft2(c) + self.dt * np.fft.fft2(Ec))
        h_hat = self._diff_if * (np.fft.fft2(h) + self.dt * np.fft.fft2(Eh))
        c_new = np.maximum(np.real(np.fft.ifft2(c_hat)), 0.0)
        h_new = np.maximum(np.real(np.fft.ifft2(h_hat)), 0.0)

        # ---- nutrient chemostat (global, shared resource) -----------------
        total_rep = float(np.sum(gross)) * self.dx**2
        dN = self.nutr_supply * (1.0 - self.Nutr) - self.nutr_use * total_rep * self.Nutr
        self.Nutr = float(np.clip(self.Nutr + self.dt * dN, 0.0, 1.0))

        self.phi, self.c, self.h = phi_new, c_new, h_new
        self._mutate()
        self.time += self.dt

    def _genotype(self):
        """g = h / c where cargo is present, clamped to the allowed range."""
        g = np.where(self.c > 1e-4, self.h / np.maximum(self.c, 1e-9),
                     0.5 * (G_MIN + G_MAX))
        return np.clip(g, G_MIN, G_MAX)

    def _mutate(self):
        """Weak mutation: diffuse the genotype and add small unbiased noise
        where cargo lives, then re-impose h = c*g (faithful-ish heredity)."""
        if self.mut <= 0:
            return
        g = self._genotype()
        # genotype spread (mutation): weak spatial smoothing + small unbiased
        # noise on live cargo only -- kept far weaker than the selection
        # differential so heredity holds and selection can act.
        g = g + 0.04 * self.mut * self._lap(g)
        live = self.c > 5e-3
        noise = 0.5 * self.mut * self.rng.standard_normal(g.shape)
        g = np.where(live, g + noise, g)
        g = np.clip(g, G_MIN, G_MAX)
        self.h = self.c * g

    def run(self, n_steps, **kw):
        for _ in range(n_steps):
            self.step(**kw)

    # -- diagnostics --------------------------------------------------------
    def mean_fitness(self, thresh=0.05):
        """Cargo-weighted population mean genotype (mean fitness)."""
        w = self.c.copy()
        w[w < thresh] = 0.0
        tot = w.sum()
        if tot <= 0:
            return 0.5 * (G_MIN + G_MAX)
        return float(np.sum(w * self._genotype()) / tot)

    def droplets(self, phi_thresh=0.0, min_area_px=12):
        """Label droplets; return list of dicts with centroid, area, fitness."""
        lab, n = label(self.phi > phi_thresh)
        g = self._genotype()
        out = []
        for i in range(1, n + 1):
            m = lab == i
            area = int(m.sum())
            if area < min_area_px:
                continue
            cy, cx = center_of_mass(m)
            cw = self.c[m].sum()
            fit = float(np.sum(self.c[m] * g[m]) / cw) if cw > 0 else np.nan
            out.append(dict(cx=cx * self.dx, cy=cy * self.dx, area=area,
                            fitness=fit, cargo=float(cw)))
        return out

    # -- fitness-coloured RGB rendering -------------------------------------
    def render_rgb(self, vmin=G_MIN, vmax=G_MAX):
        """Dense regions coloured by local fitness; dilute phase pale."""
        g = self._genotype()
        norm = Normalize(vmin, vmax)
        rgba = mpl.colormaps[FIT_CMAP](norm(g))
        dense = np.clip(0.5 * (self.phi + 1.0), 0, 1)[..., None]
        bg = np.array([0.96, 0.95, 0.92, 1.0])
        img = bg * (1 - dense) + rgba * dense
        img[..., 3] = 1.0
        return np.clip(img, 0, 1)


# ===========================================================================
# Thermal-pulse fission (chi ~ 1/T warm stripes)
# ===========================================================================

def warm_stripe_field(N, L, a0=1.0, depth=5.2, spacing=26.0, phase=0.0):
    """Cold matrix with periodic warm stripes (a dips along lines): large
    droplets spanning a stripe neck and divide; small ones survive."""
    x = np.linspace(0, L, N, endpoint=False)
    X, Y = np.meshgrid(x, x)
    band = (np.cos(2 * np.pi * (X + phase) / spacing)
            + np.cos(2 * np.pi * (Y + phase) / spacing))
    band = np.clip((band - 0.6) / 1.4, 0, 1)        # sharp warm lines
    return a0 - depth * band


# ===========================================================================
# Evolutionary driver
# ===========================================================================

def evolve(world: ProtocellWorld, generations=12, steps_per_gen=220,
           pulse_steps=90, rec_frames_per_gen=6, verbose=True):
    """Run the full Darwinian cycle; return a history dict for plotting."""
    hist = dict(gen=[], fitness=[], n_drops=[], nutrient=[], time=[],
                frames=[], frame_times=[], frame_gen=[], lineage=[],
                gen_snaps=[])
    lineage_prev = None

    def record_frame():
        hist["frames"].append(world.render_rgb())
        hist["frame_times"].append(world.time)
        hist["frame_gen"].append(len(hist["gen"]))

    # initial settle so droplets are well-formed
    record_frame()
    for g_idx in range(generations):
        # --- growth/competition phase ---
        sub = max(1, steps_per_gen // rec_frames_per_gen)
        for s in range(steps_per_gen):
            world.step(replicate=True)
            if s % sub == 0:
                record_frame()
        drops = world.droplets()
        fit = world.mean_fitness()
        hist["gen"].append(g_idx)
        hist["fitness"].append(fit)
        hist["n_drops"].append(len(drops))
        hist["nutrient"].append(world.Nutr)
        hist["time"].append(world.time)
        hist["gen_snaps"].append(world.render_rgb())

        # --- lineage bookkeeping (match daughters to nearest parent) ---
        if lineage_prev is not None:
            for d in drops:
                pi = _nearest(d, lineage_prev)
                d["parent"] = pi
        hist["lineage"].append(drops)
        lineage_prev = drops

        # --- thermal pulse: fission of the big droplets ---
        af = warm_stripe_field(world.N, world.L,
                               phase=world.rng.uniform(0, 26.0))
        for _ in range(pulse_steps):
            world.step(replicate=False, a_field=af)
            if _ % max(1, pulse_steps // 3) == 0:
                record_frame()
        # short relaxation back to the cold matrix
        world.run(40, replicate=True)
        record_frame()

        if verbose:
            print(f"  gen {g_idx:2d}:  <fitness>={fit:.3f}   "
                  f"droplets={len(drops):2d}   nutrient={world.Nutr:.2f}",
                  flush=True)
    return hist


def _nearest(d, prev):
    """Index of the nearest parent droplet by centroid distance."""
    if not prev:
        return -1
    dist = [(d["cx"] - p["cx"])**2 + (d["cy"] - p["cy"])**2 for p in prev]
    return int(np.argmin(dist))


# ===========================================================================
# Figures
# ===========================================================================

def plot_fitness_trajectory(hist, fname="darwin_fitness.png"):
    """Hero plot: mean fitness vs generation, with nutrient and droplet count,
    and the genotype distribution shifting upward."""
    set_style()
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.5, 5.6),
                                   gridspec_kw=dict(width_ratios=[1.25, 1.0]))

    gens = np.array(hist["gen"])
    fit = np.array(hist["fitness"])
    axA.plot(gens, fit, 'o-', color=ACCENT, ms=9, mec="white", mew=1.0, lw=2.6,
             zorder=6, label="population mean fitness")
    # selection-response trend (monotone rise to mutation-selection balance)
    axA.axhline(G_MAX, color=GOLD, ls=":", lw=1.6)
    axA.text(gens[-1], G_MAX - 0.04, r"$g_{\max}$ (fittest replicator)",
             ha="right", va="top", color="#8a6a1a", fontsize=10)
    axA.fill_between(gens, fit, G_MIN, color=ACCENT, alpha=0.08)
    axA.set_xlabel("generation")
    axA.set_ylabel(r"mean replication efficiency  $\langle g\rangle$")
    axA.set_ylim(G_MIN, G_MAX + 0.05)
    axA.set_title("(a) mean fitness rises by natural selection")
    axA.legend(loc="lower right", fontsize=10)

    axN = axA.twinx()
    axN.plot(gens, hist["nutrient"], 's--', color=ACCENT2, ms=5, lw=1.6,
             alpha=0.8, label="shared nutrient")
    axN.plot(gens, np.array(hist["n_drops"]) / max(hist["n_drops"]), '^:',
             color=GOLD, ms=5, lw=1.4, alpha=0.8, label="droplet count (norm.)")
    axN.set_ylabel("nutrient / droplet count", color=ACCENT2)
    axN.tick_params(axis='y', labelcolor=ACCENT2)
    axN.set_ylim(0, 1.25)
    axN.grid(False)
    axN.legend(loc="upper left", fontsize=8.5)

    # ---- (b) genotype distribution: early vs late, cargo-weighted ----------
    early = hist["lineage"][0]
    late = hist["lineage"][-1]
    eb = np.array([d["fitness"] for d in early if np.isfinite(d["fitness"])])
    ew = np.array([d["cargo"] for d in early if np.isfinite(d["fitness"])])
    lb = np.array([d["fitness"] for d in late if np.isfinite(d["fitness"])])
    lw = np.array([d["cargo"] for d in late if np.isfinite(d["fitness"])])
    bins = np.linspace(G_MIN, G_MAX, 16)
    if len(eb):
        axB.hist(eb, bins=bins, weights=ew, density=True, alpha=0.55,
                 color=ACCENT2, label=f"generation 0  ($\\langle g\\rangle$"
                                       f"={np.average(eb, weights=ew):.2f})")
    if len(lb):
        axB.hist(lb, bins=bins, weights=lw, density=True, alpha=0.6,
                 color=ACCENT, label=f"generation {gens[-1]}  "
                 f"($\\langle g\\rangle$={np.average(lb, weights=lw):.2f})")
    axB.axvline(np.average(eb, weights=ew) if len(eb) else G_MIN,
                color=ACCENT2, ls="--", lw=1.4)
    axB.axvline(np.average(lb, weights=lw) if len(lb) else G_MAX,
                color=ACCENT, ls="--", lw=1.4)
    axB.set_xlabel(r"droplet genotype  $g$ (replication efficiency)")
    axB.set_ylabel("cargo-weighted density")
    axB.set_title("(b) the population shifts toward fitter replicators")
    axB.legend(fontsize=9, loc="upper left")

    fig.suptitle("Darwinian selection among protocells: fitness increases over "
                 "generations", fontsize=14.5, y=1.01)
    fig.tight_layout()
    out = FIGURES / fname
    fig.savefig(out)
    plt.close(fig)
    print(f"  saved {out}")


def plot_generation_panels(hist, fname="darwin_generations.png", ncol=6):
    """Snapshot of the droplet field at several generations, coloured by
    fitness, showing the population becoming dominated by fit (yellow)
    droplets."""
    set_style()
    snaps = hist["gen_snaps"]
    idx = np.linspace(0, len(snaps) - 1, ncol).astype(int)
    fig, axes = plt.subplots(1, ncol, figsize=(2.7 * ncol, 3.1))
    for ax, i in zip(axes, idx):
        ax.imshow(snaps[i], origin="lower", interpolation="bilinear")
        ax.set_title(f"gen {hist['gen'][i]}\n"
                     fr"$\langle g\rangle={hist['fitness'][i]:.2f}$",
                     fontsize=10)
        ax.set_xticks([]); ax.set_yticks([])
    # colourbar for fitness
    sm = mpl.cm.ScalarMappable(Normalize(G_MIN, G_MAX), cmap=FIT_CMAP)
    cb = fig.colorbar(sm, ax=axes, fraction=0.018, pad=0.012)
    cb.set_label("droplet fitness  $g$", fontsize=10)
    fig.suptitle("Competing, growing and dividing protocells coloured by "
                 "fitness (yellow = fitter)", fontsize=13.5, y=1.06)
    out = FIGURES / fname
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out}")


def plot_lineage_tree(hist, fname="darwin_lineage.png"):
    """Lineage diagram: each droplet a node placed (generation, x-position),
    edges parent->daughter, node colour = fitness, size = cargo.  Survivors,
    deaths and divisions are visible as branching/terminating lines."""
    set_style()
    fig, ax = plt.subplots(figsize=(11.5, 6.2))
    norm = Normalize(G_MIN, G_MAX)
    cmap = mpl.colormaps[FIT_CMAP]
    lin = hist["lineage"]
    for gi, drops in enumerate(lin):
        for d in drops:
            s = 40 + 600 * d["cargo"] / (1e-6 + max(x["cargo"] for x in drops))
            ax.scatter(gi, d["cx"], s=s, color=cmap(norm(d["fitness"])),
                       ec=INK, lw=0.6, zorder=4)
            pi = d.get("parent", -1)
            if gi > 0 and pi is not None and 0 <= pi < len(lin[gi - 1]):
                p = lin[gi - 1][pi]
                ax.plot([gi - 1, gi], [p["cx"], d["cx"]], '-', color="#999",
                        lw=0.9, alpha=0.7, zorder=2)
    ax.set_xlabel("generation")
    ax.set_ylabel(r"droplet position  $x$")
    ax.set_title("Protocell lineages: survival, division and extinction\n"
                 "(node colour = fitness, size = cargo load)")
    sm = mpl.cm.ScalarMappable(norm, cmap=cmap)
    cb = fig.colorbar(sm, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("fitness  $g$")
    fig.tight_layout()
    out = FIGURES / fname
    fig.savefig(out)
    plt.close(fig)
    print(f"  saved {out}")


def animate(hist, fname="darwin_selection.mp4", fps=16):
    """Animation: droplets growing, dividing and competing, coloured by
    fitness, with a running mean-fitness readout."""
    import matplotlib.animation as animation
    set_style()
    frames = hist["frames"]
    fig, ax = plt.subplots(figsize=(5.6, 5.8))
    im = ax.imshow(frames[0], origin="lower", interpolation="bilinear")
    ax.set_xticks([]); ax.set_yticks([])
    ttl = ax.set_title("")
    sm = mpl.cm.ScalarMappable(Normalize(G_MIN, G_MAX), cmap=FIT_CMAP)
    cb = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("fitness  $g$")

    def update(f):
        im.set_data(frames[f])
        ttl.set_text(f"Darwinian protocells   $t={hist['frame_times'][f]:.0f}$"
                     f"   (gen {hist['frame_gen'][f]})")
        return im, ttl

    anim = animation.FuncAnimation(fig, update, frames=len(frames), blit=False)
    out = FIGURES / fname
    try:
        anim.save(out, writer="ffmpeg", fps=fps, dpi=120)
        print(f"  saved {out}  ({len(frames)} frames)")
    except Exception as e:
        print(f"  [skip animation: {e}]")
    plt.close(fig)


# ===========================================================================
if __name__ == "__main__":
    import sys
    quick = "--quick" in sys.argv

    print("Darwinian selection among protocells")
    print("=" * 44)
    if quick:
        world = ProtocellWorld(N=96, L=110.0, seed=1)
        world.seed_population(n_drops=7, radius=9.0)
        hist = evolve(world, generations=5, steps_per_gen=120, pulse_steps=50,
                      rec_frames_per_gen=3)
    else:
        world = ProtocellWorld(N=144, L=150.0, seed=1)
        world.seed_population(n_drops=11, radius=10.0)
        hist = evolve(world, generations=14, steps_per_gen=240, pulse_steps=95,
                      rec_frames_per_gen=6)

    print("\n[1] Fitness trajectory")
    plot_fitness_trajectory(hist)
    print("[2] Generation panels")
    plot_generation_panels(hist)
    print("[3] Lineage tree")
    plot_lineage_tree(hist)
    print("[4] Animation")
    animate(hist)
    print(f"\nDone.  fitness {hist['fitness'][0]:.3f} -> {hist['fitness'][-1]:.3f}"
          f"  over {len(hist['gen'])} generations.")
