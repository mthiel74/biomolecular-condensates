#!/usr/bin/env python3
"""Hero MP4 of 2D spinodal decomposition — clean, elegant, no clutter.

256x256 Cahn-Hilliard simulation. Just the phase field evolving.
Minimal overlay: a small time counter, nothing else. The physics
speaks for itself.
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LinearSegmentedColormap

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cahn_hilliard import CahnHilliard2D

DROP_ZONE = Path.home() / "Library/Mobile Documents/com~apple~CloudDocs/Documents/Claude"
OUT_MP4 = DROP_ZONE / "2026-06-20_condensate-spinodal-decomposition.mp4"
OUT_REPO = Path(__file__).resolve().parent.parent / "figures" / "spinodal_decomposition.mp4"

# Colormap: deep teal (dilute) → luminous gold (dense phase / condensate)
# Perceptually smooth, print-safe, more beautiful than viridis for this physics
cmap_colors = [
    (0.031, 0.075, 0.130),   # near-black teal
    (0.035, 0.145, 0.215),   # dark ocean
    (0.050, 0.230, 0.310),   # deep teal
    (0.090, 0.360, 0.390),   # teal
    (0.180, 0.500, 0.420),   # muted green-teal
    (0.380, 0.620, 0.380),   # olive transition
    (0.620, 0.720, 0.300),   # warm chartreuse
    (0.830, 0.790, 0.220),   # golden
    (0.950, 0.860, 0.310),   # bright gold
    (1.000, 0.925, 0.520),   # luminous gold
]
cmap = LinearSegmentedColormap.from_list("teal_gold", cmap_colors, N=512)

print("Spinodal decomposition animation (clean version)")
print("=" * 50)

sim = CahnHilliard2D(L=20.0, N=256, chi=1.5, N1=100, N2=1,
                     mobility=1.0, kappa=0.5, dt=0.002)
sim.initialize(phi_mean=0.3, noise_amplitude=0.01, seed=42)

total_time = 60.0
fps = 30
duration_sec = 20
n_frames = fps * duration_sec
steps_per_frame = max(1, int(total_time / sim.dt / n_frames))

print(f"  {sim.N}x{sim.N}, {n_frames} frames, {steps_per_frame} steps/frame")

frames = [sim.phi.copy()]
times = [0.0]

for i in range(n_frames):
    sim.run(steps_per_frame)
    frames.append(sim.phi.copy())
    times.append(sim.time)
    if (i + 1) % 100 == 0:
        print(f"    frame {i+1}/{n_frames}  t={sim.time:.1f}")

print("  Rendering...")

# Square figure, black background, no axes, no colorbar — just the field
px = 1080
dpi = 150
fig_size = px / dpi

fig, ax = plt.subplots(figsize=(fig_size, fig_size), facecolor='black')
fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
ax.set_position([0, 0, 1, 1])
ax.axis('off')

im = ax.imshow(frames[0], cmap=cmap, vmin=0.0, vmax=0.72,
               origin='lower', interpolation='bilinear', aspect='equal')

# Minimal time label — small, unobtrusive, bottom-right
time_text = ax.text(0.97, 0.03, "", transform=ax.transAxes,
                    ha='right', va='bottom',
                    fontsize=11, fontfamily='monospace',
                    color='white', alpha=0.55)


def update(idx):
    im.set_data(frames[idx])
    time_text.set_text(f"t = {times[idx]:.1f}")
    return [im, time_text]


ani = animation.FuncAnimation(fig, update, frames=len(frames),
                              interval=1000 / fps, blit=True)

writer = animation.FFMpegWriter(
    fps=fps, bitrate=5000,
    extra_args=['-vf', 'pad=ceil(iw/2)*2:ceil(ih/2)*2',
                '-pix_fmt', 'yuv420p'])

ani.save(str(OUT_MP4), writer=writer, dpi=dpi)
print(f"  saved {OUT_MP4}")

ani.save(str(OUT_REPO), writer=writer, dpi=dpi)
print(f"  saved {OUT_REPO}")

plt.close(fig)
print("Done.")
