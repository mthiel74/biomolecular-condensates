#!/usr/bin/env python3
"""Generate a hero MP4 animation of 2D spinodal decomposition.

Uses the CahnHilliard2D solver from cahn_hilliard.py. Produces a
smooth, labeled, publication-quality animation showing droplets
forming from a homogeneous mixture via Cahn-Hilliard dynamics.

Output: saved to iCloud drop zone as MP4.
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

# Custom colormap: deep blue (dilute) → white (critical) → warm amber (dense)
colors_list = [
    (0.05, 0.10, 0.35),  # deep navy
    (0.15, 0.30, 0.65),  # steel blue
    (0.40, 0.60, 0.85),  # soft blue
    (0.95, 0.95, 0.97),  # near-white
    (0.95, 0.75, 0.30),  # warm amber
    (0.85, 0.35, 0.10),  # burnt orange
    (0.55, 0.10, 0.05),  # deep red-brown
]
condensate_cmap = LinearSegmentedColormap.from_list("condensate", colors_list, N=256)

print("Spinodal decomposition animation")
print("=" * 50)

# Higher resolution, longer run for the hero animation
sim = CahnHilliard2D(L=20.0, N=256, chi=1.5, N1=100, N2=1,
                     mobility=1.0, kappa=0.5, dt=0.002)
sim.initialize(phi_mean=0.3, noise_amplitude=0.01, seed=42)

total_time = 60.0
fps = 30
duration_sec = 20
n_frames = fps * duration_sec
steps_per_frame = max(1, int(total_time / sim.dt / n_frames))

print(f"  Grid: {sim.N}x{sim.N}, total_time={total_time}")
print(f"  Frames: {n_frames}, steps/frame: {steps_per_frame}")
print(f"  Generating frames...")

# Collect all frames first for smooth progress
frames = []
times_arr = []
energies = []

frames.append(sim.phi.copy())
times_arr.append(sim.time)
energies.append(sim.total_free_energy())

for i in range(n_frames):
    sim.run(steps_per_frame)
    frames.append(sim.phi.copy())
    times_arr.append(sim.time)
    energies.append(sim.total_free_energy())
    if (i + 1) % 50 == 0:
        print(f"    frame {i+1}/{n_frames}  t={sim.time:.2f}")

print(f"  Simulation complete. Rendering MP4...")

# Set up the figure
fig = plt.figure(figsize=(10, 8.5), facecolor='#0a0a1a')
gs = fig.add_gridspec(5, 1, height_ratios=[0.05, 1, 0.02, 0.18, 0.02],
                      hspace=0.15, left=0.08, right=0.92, top=0.92, bottom=0.04)

ax_title = fig.add_subplot(gs[0])
ax_title.axis('off')
title_text = ax_title.text(
    0.5, 0.5,
    "Spinodal Decomposition — Biomolecular Condensate Formation",
    transform=ax_title.transAxes, ha='center', va='center',
    fontsize=16, fontweight='bold', color='white',
    fontfamily='sans-serif')

ax_main = fig.add_subplot(gs[1])
ax_energy = fig.add_subplot(gs[3])

# Main image
im = ax_main.imshow(frames[0], extent=[0, sim.L, 0, sim.L],
                    cmap=condensate_cmap, vmin=0.0, vmax=0.75,
                    origin='lower', interpolation='bilinear', aspect='equal')
ax_main.set_xlabel("x  (reduced units)", fontsize=11, color='#cccccc')
ax_main.set_ylabel("y  (reduced units)", fontsize=11, color='#cccccc')
ax_main.tick_params(colors='#888888', labelsize=9)
for spine in ax_main.spines.values():
    spine.set_color('#444444')

time_label = ax_main.text(
    0.98, 0.96, "t = 0.00", transform=ax_main.transAxes,
    ha='right', va='top', fontsize=14, fontweight='bold',
    color='white', fontfamily='monospace',
    bbox=dict(boxstyle='round,pad=0.3', facecolor='black', alpha=0.6))

# Physics annotation (appears after initial phase)
physics_label = ax_main.text(
    0.02, 0.04, "", transform=ax_main.transAxes,
    ha='left', va='bottom', fontsize=10, color='#dddddd',
    fontfamily='sans-serif',
    bbox=dict(boxstyle='round,pad=0.3', facecolor='black', alpha=0.5))

# Colorbar
cbar = fig.colorbar(im, ax=ax_main, fraction=0.03, pad=0.02)
cbar.set_label("Volume fraction  φ", fontsize=11, color='#cccccc')
cbar.ax.tick_params(colors='#888888', labelsize=9)

# Energy subplot
ax_energy.set_facecolor('#0a0a1a')
energy_line, = ax_energy.plot([], [], color='#ff9933', lw=1.5)
ax_energy.set_xlabel("Time", fontsize=9, color='#aaaaaa')
ax_energy.set_ylabel("Free energy F[φ]", fontsize=9, color='#aaaaaa')
ax_energy.tick_params(colors='#666666', labelsize=8)
for spine in ax_energy.spines.values():
    spine.set_color('#333333')
ax_energy.set_xlim(0, total_time)
e_arr = np.array(energies)
ax_energy.set_ylim(e_arr.min() * 1.05, e_arr.max() * 0.98)
ax_energy.grid(True, alpha=0.15, color='#444444')

def get_phase_text(t):
    if t < 1:
        return "Homogeneous mixture — metastable"
    elif t < 5:
        return "Spinodal instability — uphill diffusion"
    elif t < 15:
        return "Droplet formation — phase separation"
    elif t < 35:
        return "Coarsening — Ostwald ripening (⟨R⟩ ~ t¹ᐟ³)"
    else:
        return "Late-stage coarsening — fewer, larger droplets"


def update(frame_idx):
    im.set_data(frames[frame_idx])
    time_label.set_text(f"t = {times_arr[frame_idx]:.2f}")
    physics_label.set_text(get_phase_text(times_arr[frame_idx]))
    energy_line.set_data(times_arr[:frame_idx+1], energies[:frame_idx+1])
    return [im, time_label, physics_label, energy_line]


ani = animation.FuncAnimation(fig, update, frames=len(frames),
                              interval=1000/fps, blit=True)

writer = animation.FFMpegWriter(fps=fps, bitrate=4000,
                                extra_args=['-vf', 'pad=ceil(iw/2)*2:ceil(ih/2)*2',
                                            '-pix_fmt', 'yuv420p'])
ani.save(str(OUT_MP4), writer=writer, dpi=150)
print(f"  saved {OUT_MP4}")

ani.save(str(OUT_REPO), writer=writer, dpi=150)
print(f"  saved {OUT_REPO}")

plt.close(fig)
print("Done.")
