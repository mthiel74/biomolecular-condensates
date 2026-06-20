#!/usr/bin/env python3
"""Animate protocell population dynamics — proto-Darwinian selection.

Shows a population of condensate droplets growing, dividing, and
competing. Higher-K droplets (warmer colors) out-compete lower-K
ones (cooler colors). Tracks mean fitness over time.

Output: saved to iCloud drop zone as MP4.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from pathlib import Path

DROP_ZONE = Path.home() / "Library/Mobile Documents/com~apple~CloudDocs/Documents/Claude"
OUT_MP4 = DROP_ZONE / "2026-06-20_condensate-protocell-selection.mp4"
OUT_REPO = Path(__file__).resolve().parent.parent / "figures" / "protocell_selection.mp4"

print("Protocell selection animation")
print("=" * 50)


class Protocell:
    def __init__(self, x, y, radius, K, gen=0):
        self.x = x
        self.y = y
        self.radius = radius
        self.K = K
        self.gen = gen
        self.alive = True

    def growth_rate(self):
        return 0.002 * self.K

    def grow(self, dt, resource_factor=1.0):
        self.radius += self.growth_rate() * dt * resource_factor
        self.radius = max(self.radius, 0.02)


def simulate_protocells(n_init=30, n_steps=800, box_size=10.0,
                        r_divide=0.6, r_death=0.08, max_pop=80, seed=42):
    """Run protocell simulation, return frame data for animation."""
    rng = np.random.default_rng(seed)

    cells = []
    for _ in range(n_init):
        K = rng.uniform(10, 80)
        r = rng.uniform(0.15, 0.35)
        x = rng.uniform(r, box_size - r)
        y = rng.uniform(r, box_size - r)
        cells.append(Protocell(x, y, r, K))

    frame_data = []

    for step in range(n_steps):
        resource = max(0.3, 1.0 - len(cells) / max_pop)

        for c in cells:
            c.grow(1.0, resource)

        # Division
        new_cells = []
        for c in list(cells):
            if c.radius >= r_divide and len(cells) + len(new_cells) < max_pop:
                d_angle = rng.uniform(0, 2 * np.pi)
                offset = c.radius * 0.4
                r_new = c.radius * 0.55
                K_mut = c.K + rng.normal(0, 2.0)
                K_mut = max(5, K_mut)

                c.radius = r_new
                daughter = Protocell(
                    np.clip(c.x + offset * np.cos(d_angle), 0.5, box_size - 0.5),
                    np.clip(c.y + offset * np.sin(d_angle), 0.5, box_size - 0.5),
                    r_new, K_mut, c.gen + 1)
                c.gen += 1
                new_cells.append(daughter)
        cells.extend(new_cells)

        # Death — low-K cells more likely to die under competition
        for c in cells:
            if c.radius < r_death:
                c.alive = False
            elif resource < 0.5 and rng.random() < 0.01 * (80 / max(c.K, 1)):
                c.alive = False

        cells = [c for c in cells if c.alive]

        # Slight random drift
        for c in cells:
            c.x = np.clip(c.x + rng.normal(0, 0.03), c.radius, box_size - c.radius)
            c.y = np.clip(c.y + rng.normal(0, 0.03), c.radius, box_size - c.radius)

        if step % 2 == 0:
            frame_data.append({
                'step': step,
                'xs': [c.x for c in cells],
                'ys': [c.y for c in cells],
                'rs': [c.radius for c in cells],
                'Ks': [c.K for c in cells],
                'gens': [c.gen for c in cells],
                'n': len(cells),
                'mean_K': np.mean([c.K for c in cells]) if cells else 0,
                'max_gen': max((c.gen for c in cells), default=0),
            })

    return frame_data


print("  Simulating protocell population...")
frame_data = simulate_protocells()
print(f"  {len(frame_data)} frames captured")

# Build the animation
fig = plt.figure(figsize=(12, 7), facecolor='#0a0e1a')
gs = fig.add_gridspec(1, 2, width_ratios=[1.3, 1], wspace=0.25,
                      left=0.06, right=0.96, top=0.88, bottom=0.10)

ax_field = fig.add_subplot(gs[0])
ax_stats = fig.add_subplot(gs[1])

fig.suptitle("Proto-Darwinian Selection in Condensate Populations",
             fontsize=15, fontweight='bold', color='white', y=0.95)

# Field setup
ax_field.set_xlim(0, 10)
ax_field.set_ylim(0, 10)
ax_field.set_aspect('equal')
ax_field.set_facecolor('#060a18')
ax_field.set_xlabel("x", color='#999999', fontsize=10)
ax_field.set_ylabel("y", color='#999999', fontsize=10)
ax_field.tick_params(colors='#666666', labelsize=8)
for spine in ax_field.spines.values():
    spine.set_color('#333355')

# K color normalization
K_min, K_max = 5, 100

# Stats panel setup
ax_stats.set_facecolor('#0a0e1a')
for spine in ax_stats.spines.values():
    spine.set_color('#333355')
ax_stats.tick_params(colors='#666666', labelsize=8)

# Pre-compute stats arrays
steps_arr = [f['step'] for f in frame_data]
meanK_arr = [f['mean_K'] for f in frame_data]
n_arr = [f['n'] for f in frame_data]
gen_arr = [f['max_gen'] for f in frame_data]

# Stats lines
ax_K = ax_stats
ax_K.set_xlim(0, steps_arr[-1])
ax_K.set_ylim(0, 90)
ax_K.set_xlabel("Time step", color='#999999', fontsize=10)
ax_K.set_ylabel("Mean partition coefficient K", color='#ff9944', fontsize=10)
line_K, = ax_K.plot([], [], color='#ff9944', lw=2, label='Mean K')

ax_n = ax_K.twinx()
ax_n.set_ylim(0, 90)
ax_n.set_ylabel("Population size", color='#44aaff', fontsize=10)
ax_n.tick_params(colors='#666666', labelsize=8)
line_n, = ax_n.plot([], [], color='#44aaff', lw=1.5, alpha=0.7, label='Population')
ax_n.spines['right'].set_color('#333355')

ax_K.legend(handles=[line_K, line_n], loc='upper left', fontsize=9,
            facecolor='#0a0e1a', edgecolor='#333355', labelcolor='#cccccc')

step_label = ax_field.text(
    0.98, 0.97, "", transform=ax_field.transAxes,
    ha='right', va='top', fontsize=12, fontweight='bold',
    color='white', fontfamily='monospace',
    bbox=dict(boxstyle='round,pad=0.3', facecolor='black', alpha=0.6))

pop_label = ax_field.text(
    0.02, 0.97, "", transform=ax_field.transAxes,
    ha='left', va='top', fontsize=10,
    color='#aaccff', fontfamily='sans-serif',
    bbox=dict(boxstyle='round,pad=0.3', facecolor='black', alpha=0.5))


def update(idx):
    f = frame_data[idx]
    for p in ax_field.patches[:]:
        p.remove()

    for x, y, r, K in zip(f['xs'], f['ys'], f['rs'], f['Ks']):
        K_norm = np.clip((K - K_min) / (K_max - K_min), 0, 1)
        # Cool blue (low K) → warm orange (high K)
        color = plt.cm.plasma(K_norm)
        alpha = 0.5 + 0.4 * K_norm
        circle = Circle((x, y), r, facecolor=(*color[:3], alpha),
                        edgecolor=(*color[:3], 0.9), linewidth=0.8)
        ax_field.add_patch(circle)

    step_label.set_text(f"t = {f['step']}")
    pop_label.set_text(f"n = {f['n']}   ⟨K⟩ = {f['mean_K']:.1f}")

    line_K.set_data(steps_arr[:idx+1], meanK_arr[:idx+1])
    line_n.set_data(steps_arr[:idx+1], n_arr[:idx+1])

    return []


fps = 30
print(f"  Rendering {len(frame_data)} frames at {fps} fps...")

ani = animation.FuncAnimation(fig, update, frames=len(frame_data),
                              interval=1000/fps, blit=False)

writer = animation.FFMpegWriter(fps=fps, bitrate=3000,
                                extra_args=['-vf', 'pad=ceil(iw/2)*2:ceil(ih/2)*2',
                                            '-pix_fmt', 'yuv420p'])
ani.save(str(OUT_MP4), writer=writer, dpi=120)
print(f"  saved {OUT_MP4}")

ani.save(str(OUT_REPO), writer=writer, dpi=120)
print(f"  saved {OUT_REPO}")

plt.close(fig)
print("Done.")
