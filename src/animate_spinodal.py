#!/usr/bin/env python3
"""Hero MP4 of 2D spinodal decomposition — clean, elegant, no clutter.

Uses a polynomial double-well potential f(φ)=φ²(1-φ)² for numerical
stability and visual beauty. Semi-implicit spectral solver with exact
mass conservation (no clipping artefacts). 256x256, 30fps, 15 seconds.

The FH-based solver in cahn_hilliard.py is used for the static physics
figures; this animation uses the simpler potential to produce the most
visually dramatic phase separation.
"""

from pathlib import Path
import numpy as np
from numpy.fft import fft2, ifft2, fftfreq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LinearSegmentedColormap

DROP_ZONE = Path.home() / "Library/Mobile Documents/com~apple~CloudDocs/Documents/Claude"
OUT_MP4 = DROP_ZONE / "2026-06-20_condensate-spinodal-decomposition.mp4"
OUT_REPO = Path(__file__).resolve().parent.parent / "figures" / "spinodal_decomposition.mp4"


def mu_doublewell(phi):
    """Chemical potential from f(φ)=φ²(1-φ)², μ=2φ(1-φ)(1-2φ)."""
    return 2.0 * phi * (1.0 - phi) * (1.0 - 2.0 * phi)


def run_simulation(N=256, L=20.0, kappa=0.5, M=1.0, dt=0.005,
                   phi_mean=0.35, noise=0.05, total_time=50.0,
                   frames=450, seed=42):
    """Run Cahn-Hilliard with polynomial double-well. Returns frame list."""
    dx = L / N
    kx = 2 * np.pi * fftfreq(N, d=dx)
    KX, KY = np.meshgrid(kx, kx)
    k2 = KX**2 + KY**2

    # Stabilisation: A ≥ max|d²f/dφ²|. For f=φ²(1-φ)²,
    # d²f/dφ² = 12φ²-12φ+2, max at φ=0 or 1: value=2. Use A=3.
    A = 3.0
    denom = 1.0 + dt * M * k2 * (A + kappa * k2)

    rng = np.random.default_rng(seed)
    phi = phi_mean + noise * rng.standard_normal((N, N))
    phi = np.clip(phi, 0.01, 0.99)

    steps_per_frame = max(1, int(total_time / dt / frames))
    print(f"  {N}x{N}, {frames} frames, {steps_per_frame} steps/frame, dt={dt}")

    snapshots = [phi.copy()]
    times = [0.0]
    t = 0.0

    for f_idx in range(frames):
        for _ in range(steps_per_frame):
            mu = mu_doublewell(phi)
            g = mu - A * phi
            g_hat = fft2(g)
            phi_hat = fft2(phi)
            phi_hat = (phi_hat - dt * M * k2 * g_hat) / denom
            phi = np.real(ifft2(phi_hat))
            t += dt

        snapshots.append(phi.copy())
        times.append(t)
        if (f_idx + 1) % 50 == 0:
            print(f"    frame {f_idx+1}/{frames}  t={t:.1f}"
                  f"  min={phi.min():.3f} max={phi.max():.3f}"
                  f"  mean={phi.mean():.4f}")

    return snapshots, times


# Colormap: dark indigo (dilute) → cream/gold (dense condensate)
cmap_colors = [
    (0.02, 0.02, 0.10),    # near-black indigo
    (0.06, 0.05, 0.22),    # deep indigo
    (0.12, 0.10, 0.38),    # indigo
    (0.22, 0.15, 0.52),    # purple-indigo
    (0.35, 0.22, 0.60),    # warm purple
    (0.52, 0.32, 0.58),    # mauve
    (0.72, 0.48, 0.45),    # dusty rose
    (0.88, 0.68, 0.38),    # warm sand
    (0.96, 0.85, 0.55),    # cream-gold
    (1.00, 0.96, 0.80),    # bright cream
]
cmap = LinearSegmentedColormap.from_list("indigo_gold", cmap_colors, N=512)


if __name__ == "__main__":
    print("Spinodal decomposition animation")
    print("=" * 50)

    snapshots, times = run_simulation(
        N=256, L=20.0, kappa=0.5, M=1.0, dt=0.005,
        phi_mean=0.35, noise=0.05, total_time=50.0,
        frames=450, seed=42)

    print(f"  Mass conservation check: initial mean={snapshots[0].mean():.6f}"
          f"  final mean={snapshots[-1].mean():.6f}")

    # Square figure, black background, no axes — just the field
    px = 1080
    dpi = 150
    fig_size = px / dpi

    fig, ax = plt.subplots(figsize=(fig_size, fig_size), facecolor='black')
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    ax.set_position([0, 0, 1, 1])
    ax.axis('off')

    im = ax.imshow(snapshots[0], cmap=cmap, vmin=0.0, vmax=1.0,
                   origin='lower', interpolation='bilinear', aspect='equal')

    time_text = ax.text(0.97, 0.03, "", transform=ax.transAxes,
                        ha='right', va='bottom',
                        fontsize=11, fontfamily='monospace',
                        color='white', alpha=0.45)

    fps = 30

    def update(idx):
        im.set_data(snapshots[idx])
        time_text.set_text(f"t = {times[idx]:.1f}")
        return [im, time_text]

    ani = animation.FuncAnimation(fig, update, frames=len(snapshots),
                                  interval=1000 / fps, blit=True)

    writer = animation.FFMpegWriter(
        fps=fps, bitrate=5000,
        extra_args=['-vf', 'pad=ceil(iw/2)*2:ceil(ih/2)*2',
                    '-pix_fmt', 'yuv420p'])

    print("  Rendering MP4...")
    ani.save(str(OUT_MP4), writer=writer, dpi=dpi)
    print(f"  saved {OUT_MP4}")
    ani.save(str(OUT_REPO), writer=writer, dpi=dpi)
    print(f"  saved {OUT_REPO}")

    plt.close(fig)

    # Verify: extract 3 frames with ffmpeg and save as PNGs
    import subprocess
    for t_sec, label in [(0, "start"), (7, "mid"), (14, "end")]:
        out_png = f"/tmp/verify_{label}.png"
        subprocess.run(["ffmpeg", "-y", "-ss", str(t_sec), "-i", str(OUT_REPO),
                        "-frames:v", "1", out_png],
                       capture_output=True)
    print("  Verification frames saved to /tmp/verify_*.png")
    print("Done.")
