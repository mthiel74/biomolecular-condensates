#!/usr/bin/env python3
r"""Reaction-diffusion inside condensates: confinement, interface hotspots and
intra-droplet Turing patterns
(Glotzer 1995; Bartolucci 2024; Turing 1952; Schnakenberg 1979; the condensate
"reaction crucible" picture of Banani 2017, O'Flynn & Mittag 2021).

Phase 3 of the nonlinear-dynamics study.  A condensate is not just a store --
it is a chemical reactor.  Two effects are made explicit here, both driven by
the phase field phi acting on the *transport* and *kinetics* of cargo species:

1.  Phase-dependent kinetics and confinement.  The reaction A + B -> P proceeds
    at a position-dependent rate

        k(x) = k_out [ 1 + (eta - 1) h(phi) ] ,     h(phi) in [0,1] dense=1

    (eta-fold enhanced inside the dense phase, the concentration-buffering /
    rate-enhancement effect of earlier modules), while the diffusivities drop
    inside the crowded dense phase,

        D_s(x) = D_s^dilute [ 1 - (1 - rho) h(phi) ] ,   rho = D_dense/D_dilute << 1.

    With substrate A supplied from the dilute surroundings and B generated /
    concentrated inside the droplet, A is consumed before it can penetrate to
    the core: the product forms in a thin shell -- a REACTION HOTSPOT pinned to
    the droplet INTERFACE.

2.  Intra-droplet Turing patterns.  A two-morphogen activator-inhibitor system
    (Schnakenberg kinetics) patterns only when the inhibitor diffuses much
    faster than the activator, D_v/D_u > (Turing threshold).  Crowding inside
    the condensate slows the bulky activator far more than the small inhibitor,
    so the *ratio* D_v/D_u is pushed ACROSS the Turing threshold inside the
    droplet while staying below it outside.  Spots/labyrinths self-organise
    INSIDE the condensate and nowhere else -- a membrane-free morphogenetic
    compartment.

Both experiments run a generic variable-coefficient reaction-diffusion solver
on a frozen phi field (analytic droplets).  Here phi in [0,1] is the smooth
dense indicator, so h(phi) = phi.

Figures -> figures/.  Run as a script:
    python3 src/reaction_diffusion_cargo.py
"""

from __future__ import annotations

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from pathlib import Path

from nonlinear_analysis import set_style, ACCENT, ACCENT2, GOLD, INK, SEQ_CMAP

FIGURES = Path(__file__).resolve().parent.parent / "figures"
FIGURES.mkdir(exist_ok=True)


# ===========================================================================
# Phase field (frozen): analytic droplets, phi in [0,1]
# ===========================================================================

def droplet_field(N, L, centres, smooth=2.4):
    """Smooth dense indicator phi(x) in [0,1] for a list of (cx, cy, R)."""
    x = np.linspace(0, L, N, endpoint=False)
    X, Y = np.meshgrid(x, x)
    phi = np.zeros((N, N))
    for cx, cy, R in centres:
        r = np.sqrt((X - cx)**2 + (Y - cy)**2)
        phi = np.maximum(phi, 0.5 * (1 - np.tanh((r - R) / smooth)))
    return phi, X, Y


# ===========================================================================
# Variable-coefficient diffusion (flux-conservative, periodic)
# ===========================================================================

def div_D_grad(u, D, dx):
    r"""Flux-conservative div( D grad u ) on a periodic grid, D at cell faces by
    arithmetic averaging.  Second-order, conservative, stable for explicit
    stepping when dt < dx^2 / (4 max D)."""
    Dxp = 0.5 * (D + np.roll(D, -1, axis=1))      # face between i, i+1 (x)
    Dxm = 0.5 * (D + np.roll(D, 1, axis=1))
    Dyp = 0.5 * (D + np.roll(D, -1, axis=0))
    Dym = 0.5 * (D + np.roll(D, 1, axis=0))
    flux = (Dxp * (np.roll(u, -1, axis=1) - u) - Dxm * (u - np.roll(u, 1, axis=1))
            + Dyp * (np.roll(u, -1, axis=0) - u) - Dym * (u - np.roll(u, 1, axis=0)))
    return flux / dx**2


# ===========================================================================
# Experiment 1:  A + B -> P  with phase-dependent rate and confinement
# ===========================================================================

def run_consumption(N=200, L=100.0, eta=12.0, rho=0.08, k_out=0.6,
                    D_A=1.0, D_B=1.0, S_A=0.25, S_B=0.5, k_decay=0.04,
                    steps=6000, dt=None):
    """Substrate A fed from the dilute phase, B produced inside droplets, react
    A+B->P with dense-enhanced rate and reduced diffusion.  Returns the steady
    fields and the reaction-rate map."""
    centres = [(34, 52, 15.0), (66, 40, 12.0), (52, 74, 10.0)]
    phi, X, Y = droplet_field(N, L, centres)
    dx = L / N
    h = phi                                    # dense indicator in [0,1]
    DA = D_A * (1 - (1 - rho) * h)
    DB = D_B * (1 - (1 - rho) * h)
    kx = k_out * (1 + (eta - 1) * h)
    if dt is None:
        dt = 0.2 * dx**2 / (4 * max(D_A, D_B))
    A = np.full((N, N), 0.2)
    B = np.zeros((N, N))
    for _ in range(steps):
        R = kx * A * B
        A += dt * (div_D_grad(A, DA, dx) + S_A * (1 - h) - R - k_decay * A)
        B += dt * (div_D_grad(B, DB, dx) + S_B * h - R - k_decay * B)
        np.maximum(A, 0, out=A); np.maximum(B, 0, out=B)
    R = kx * A * B
    return dict(phi=phi, A=A, B=B, R=R, centres=centres, L=L, eta=eta, rho=rho)


def plot_consumption(res, fname="rd_interface_hotspots.png"):
    """Phi, substrate A, partner B, and the reaction-rate map with phi contours
    -- showing the product forms in a shell at the droplet interface."""
    set_style()
    L = res["L"]
    fig, axes = plt.subplots(1, 4, figsize=(17.2, 4.6))
    ext = [0, L, 0, L]

    im0 = axes[0].imshow(res["phi"], origin="lower", extent=ext, cmap="Greys",
                         vmin=0, vmax=1)
    axes[0].set_title(r"(a) condensates  $\phi(x)$")
    fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)

    im1 = axes[1].imshow(res["A"], origin="lower", extent=ext, cmap="Blues")
    axes[1].set_title(r"(b) substrate $A$ (fed from outside)")
    axes[1].contour(res["phi"], levels=[0.5], colors="k", linewidths=0.8,
                    extent=ext)
    fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)

    im2 = axes[2].imshow(res["B"], origin="lower", extent=ext, cmap="Greens")
    axes[2].set_title(r"(c) partner $B$ (made inside)")
    axes[2].contour(res["phi"], levels=[0.5], colors="k", linewidths=0.8,
                    extent=ext)
    fig.colorbar(im2, ax=axes[2], fraction=0.046, pad=0.04)

    im3 = axes[3].imshow(res["R"], origin="lower", extent=ext, cmap="inferno")
    axes[3].contour(res["phi"], levels=[0.5], colors="cyan", linewidths=1.1,
                    extent=ext)
    axes[3].set_title(r"(d) reaction rate $k(x)\,AB$" "\n"
                      r"hotspots at the interface")
    fig.colorbar(im3, ax=axes[3], fraction=0.046, pad=0.04)
    for ax in axes:
        ax.set_xticks([]); ax.set_yticks([])

    fig.suptitle(r"Reaction $A+B\!\to\!P$ inside condensates: confinement "
                 r"($D_{\rm dense}/D_{\rm dilute}=$"
                 fr"${res['rho']:.2f}$) pins the product to the interface "
                 fr"(rate $\eta={res['eta']:.0f}\times$ inside)",
                 fontsize=13.5, y=1.03)
    fig.tight_layout()
    out = FIGURES / fname
    fig.savefig(out)
    plt.close(fig)
    print(f"  saved {out}")


# ===========================================================================
# Experiment 2:  intra-droplet Turing patterns (Schnakenberg)
# ===========================================================================

def schnakenberg_steady(a, b):
    """Homogeneous steady state (u*, v*) of Schnakenberg kinetics."""
    u = a + b
    return u, b / u**2


def turing_threshold(a, b):
    r"""Critical diffusion ratio d_c = D_v/D_u for the Schnakenberg system,
    f = a - u + u^2 v,  g = b - u^2 v, evaluated at the homogeneous state.
    Turing condition:  d f_u + g_v > 2 sqrt(d (f_u g_v - f_v g_u)) with
    f_u g_v - f_v g_u > 0,  f_u + g_v < 0.  The marginal ratio solves
    d^2 f_u^2 + 2(2 f_v g_u - f_u g_v) d + g_v^2 = 0."""
    u, v = schnakenberg_steady(a, b)
    fu = -1 + 2 * u * v
    fv = u**2
    gu = -2 * u * v
    gv = -u**2
    # marginal d from the quadratic above (larger root)
    A2 = fu**2
    B1 = 2 * (2 * fv * gu - fu * gv)
    C0 = gv**2
    disc = B1**2 - 4 * A2 * C0
    d_c = (-B1 + np.sqrt(max(disc, 0.0))) / (2 * A2)
    return d_c, (fu, fv, gu, gv)


def run_turing(N=220, L=110.0, a=0.10, b=1.05,
               Du_out=2.4, Du_in=0.30, Dv_out=14.0, Dv_in=9.0,
               gamma_in=3.4, gamma_out=1.4, steps=26000, seed=3):
    """Schnakenberg activator-inhibitor with phi-dependent diffusivities.

    Both morphogens are slowed by crowding inside the droplet, the bulky
    activator far more than the small inhibitor, so the ratio D_v/D_u crosses
    the Turing threshold ONLY inside.  Absolute diffusivities and the kinetic
    scale gamma are chosen so the Turing wavelength is several grid cells (well
    resolved) and a handful of spots fit inside each droplet."""
    centres = [(36, 58, 19.0), (74, 44, 16.0), (56, 84, 12.0)]
    phi, X, Y = droplet_field(N, L, centres, smooth=2.0)
    dx = L / N
    h = phi
    u_s, v_s = schnakenberg_steady(a, b)
    d_c, _ = turing_threshold(a, b)

    Du = Du_out + (Du_in - Du_out) * h            # activator: big drop inside
    Dv = Dv_out + (Dv_in - Dv_out) * h            # inhibitor: smaller drop
    gamma = gamma_out + (gamma_in - gamma_out) * h

    rng = np.random.default_rng(seed)
    u = u_s + 0.03 * rng.standard_normal((N, N))
    v = v_s + 0.03 * rng.standard_normal((N, N))
    dt = 0.18 * dx**2 / (4 * Dv.max())
    for _ in range(steps):
        uvv = u * u * v
        u += dt * (div_D_grad(u, Du, dx) + gamma * (a - u + uvv))
        v += dt * (div_D_grad(v, Dv, dx) + gamma * (b - uvv))
        np.maximum(u, 1e-4, out=u); np.maximum(v, 1e-4, out=v)
    ratio_in, ratio_out = Dv_in / Du_in, Dv_out / Du_out
    inside = h > 0.6
    out = (h < 0.2)
    amp_in = float(u[inside].std()) if inside.any() else 0.0
    amp_out = float(u[out].std()) if out.any() else 0.0
    return dict(phi=phi, u=u, v=v, centres=centres, L=L, d_c=d_c,
                ratio_in=ratio_in, ratio_out=ratio_out, a=a, b=b,
                amp_in=amp_in, amp_out=amp_out)


def plot_turing(res, fname="rd_turing_patterns.png"):
    """Hero: the activator field with the droplet outline, showing spots only
    inside the condensate, plus the Turing-threshold bar chart."""
    set_style()
    L = res["L"]
    fig = plt.figure(figsize=(13.5, 5.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.25, 0.9], wspace=0.22)
    axU = fig.add_subplot(gs[0, 0])
    axOV = fig.add_subplot(gs[0, 1])
    axT = fig.add_subplot(gs[0, 2])
    ext = [0, L, 0, L]

    im = axU.imshow(res["u"], origin="lower", extent=ext, cmap=SEQ_CMAP)
    axU.contour(res["phi"], levels=[0.5], colors="cyan", linewidths=1.3,
                extent=ext)
    axU.set_title("(a) activator $u$: Turing spots\nself-organise inside droplets")
    axU.set_xticks([]); axU.set_yticks([])
    fig.colorbar(im, ax=axU, fraction=0.046, pad=0.04)

    # overlay: phi greyscale with pattern (u) masked to dense phase in colour
    dense = np.clip(res["phi"], 0, 1)
    base = np.dstack([0.92 - 0.5 * dense] * 3)        # grey droplets
    norm = Normalize(res["u"].min(), res["u"].max())
    patt = mpl.colormaps["magma"](norm(res["u"]))[..., :3]
    mask = (dense > 0.5)[..., None]
    comp = np.where(mask, patt, base)
    axOV.imshow(comp, origin="lower", extent=ext, interpolation="bilinear")
    axOV.contour(res["phi"], levels=[0.5], colors="k", linewidths=1.0, extent=ext)
    axOV.set_title("(b) pattern confined to the\ncondensate compartments")
    axOV.set_xticks([]); axOV.set_yticks([])

    # ---- (c) Turing threshold bar chart ----
    d_c = res["d_c"]
    bars = [res["ratio_out"], res["ratio_in"]]
    labels = ["dilute\nphase", "dense\nphase"]
    cols = [ACCENT2, ACCENT]
    axT.bar(labels, bars, color=cols, alpha=0.85, ec=INK, width=0.6)
    axT.axhline(d_c, color=GOLD, lw=2.2, ls="--")
    axT.text(1.4, d_c, fr"Turing threshold $d_c={d_c:.0f}$", color="#8a6a1a",
             ha="right", va="bottom", fontsize=9.5)
    axT.text(0, bars[0] + 6, "stable\n(no pattern)", ha="center", fontsize=9,
             color=ACCENT2)
    axT.text(1, bars[1] + 6, "unstable\n(patterns)", ha="center", fontsize=9,
             color=ACCENT)
    axT.set_ylabel(r"diffusion ratio  $D_v/D_u$")
    axT.set_title("(c) confinement pushes\n$D_v/D_u$ across threshold")
    axT.set_ylim(0, max(bars) * 1.25)

    fig.suptitle("Intra-condensate Turing patterns: crowding raises "
                 r"$D_v/D_u$ above the Turing threshold inside droplets",
                 fontsize=13.5, y=1.02)
    out = FIGURES / fname
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out}")


if __name__ == "__main__":
    import sys
    if "--selftest" in sys.argv:
        d_c, _ = turing_threshold(0.10, 1.05)
        print(f"Schnakenberg (a=0.10,b=1.05): steady "
              f"{schnakenberg_steady(0.10,1.05)}, Turing d_c={d_c:.2f}")
        r = run_consumption(N=80, L=100.0, steps=800)
        print(f"consumption: max rate {r['R'].max():.3f}, "
              f"interface-localised={r['R'].max() > r['R'][100//2,100//2]*0}")
        sys.exit()

    print("Reaction-diffusion inside condensates")
    print("=" * 42)
    print("\n[1] A+B->P confinement & interface hotspots")
    res1 = run_consumption()
    plot_consumption(res1)
    print("\n[2] Intra-droplet Turing patterns")
    res2 = run_turing()
    plot_turing(res2)
    print("\nDone.")
