#!/usr/bin/env python3
"""Animate protocell population dynamics — proto-Darwinian selection.

Shows condensate droplets growing, dividing, and competing.
Higher-K droplets (warm colors) out-compete lower-K ones (cool).
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pathlib import Path

DROP_ZONE = Path.home() / "Library/Mobile Documents/com~apple~CloudDocs/Documents/Claude"
OUT_MP4 = DROP_ZONE / "2026-06-20_condensate-protocell-selection.mp4"
OUT_REPO = Path(__file__).resolve().parent.parent / "figures" / "protocell_selection.mp4"


def simulate(n_init=30, n_steps=1200, box=10.0,
             r_divide=0.55, r_death=0.06, max_pop=80, seed=42):
    rng = np.random.default_rng(seed)
    xs = rng.uniform(0.8, box - 0.8, n_init)
    ys = rng.uniform(0.8, box - 0.8, n_init)
    rs = rng.uniform(0.15, 0.30, n_init)
    Ks = rng.uniform(10, 80, n_init)

    frames = []
    for step in range(n_steps):
        n = len(xs)
        if n == 0:
            break
        resource = max(0.3, 1.0 - n / max_pop)

        # Grow — rate proportional to K, slow enough for visible dynamics
        growth = 0.0004 * Ks * resource
        rs = rs + growth
        rs = np.maximum(rs, 0.02)

        # Division
        dividers = (rs >= r_divide) & (n + 1 <= max_pop)
        n_div = dividers.sum()
        if n_div > 0:
            # Limit divisions per step to avoid population explosion
            if n + n_div > max_pop:
                idxs = np.where(dividers)[0]
                rng.shuffle(idxs)
                keep = max_pop - n
                dividers[:] = False
                dividers[idxs[:keep]] = True
                n_div = keep

            ang = rng.uniform(0, 2 * np.pi, n_div)
            off = rs[dividers] * 0.5
            r_new = rs[dividers] * 0.52
            K_new = np.maximum(5, Ks[dividers] + rng.normal(0, 2.0, n_div))

            dx = off * np.cos(ang)
            dy = off * np.sin(ang)

            new_xs = np.clip(xs[dividers] + dx, 0.5, box - 0.5)
            new_ys = np.clip(ys[dividers] + dy, 0.5, box - 0.5)

            rs[dividers] = r_new  # parent shrinks
            xs = np.concatenate([xs, new_xs])
            ys = np.concatenate([ys, new_ys])
            rs = np.concatenate([rs, r_new.copy()])
            Ks = np.concatenate([Ks, K_new])

        # Death
        alive = rs >= r_death
        # Competition pressure on low-K cells
        if resource < 0.5:
            death_prob = 0.008 * (80.0 / np.maximum(Ks, 1.0))
            alive &= rng.random(len(xs)) > death_prob

        xs, ys, rs, Ks = xs[alive], ys[alive], rs[alive], Ks[alive]

        # Drift
        xs = np.clip(xs + rng.normal(0, 0.02, len(xs)), rs, box - rs)
        ys = np.clip(ys + rng.normal(0, 0.02, len(ys)), rs, box - rs)

        if step % 3 == 0:
            frames.append({
                'step': step,
                'xs': xs.copy(), 'ys': ys.copy(),
                'rs': rs.copy(), 'Ks': Ks.copy(),
                'n': len(xs),
                'mean_K': float(np.mean(Ks)) if len(Ks) > 0 else 0,
            })
    return frames


print("Protocell selection animation")
print("=" * 50)
print("  Simulating...")
frames = simulate()
print(f"  {len(frames)} frames, final pop={frames[-1]['n']}, "
      f"final ⟨K⟩={frames[-1]['mean_K']:.1f}")

# Figure
fig = plt.figure(figsize=(12, 6.5), facecolor='#080c18')
gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1], wspace=0.22,
                      left=0.05, right=0.97, top=0.88, bottom=0.12)

ax_f = fig.add_subplot(gs[0])
ax_s = fig.add_subplot(gs[1])

fig.suptitle("Proto-Darwinian Selection in Condensate Populations",
             fontsize=14, fontweight='bold', color='white', y=0.95)

ax_f.set_xlim(0, 10)
ax_f.set_ylim(0, 10)
ax_f.set_aspect('equal')
ax_f.set_facecolor('#060a18')
ax_f.tick_params(colors='#555555', labelsize=8)
for s in ax_f.spines.values():
    s.set_color('#222244')

K_lo, K_hi = 5, 100

# Scatter — convert radius in data units to marker size in points²
# Figure is 12in wide, left panel ~57% → ~6.84in, data range 10 → 0.684 in/unit
# At dpi=120 → 82 px/unit. radius 0.3 → 24.6px → area ~ 1900 pt²
def r_to_s(r):
    px_per_unit = 82.0
    return (np.minimum(r, 0.5) * px_per_unit) ** 2 * np.pi / 4

f0 = frames[0]
K_n = np.clip((f0['Ks'] - K_lo) / (K_hi - K_lo), 0, 1)
scatter = ax_f.scatter(f0['xs'], f0['ys'], s=r_to_s(f0['rs']),
                       c=K_n, cmap='plasma', vmin=0, vmax=1,
                       alpha=0.8, edgecolors='white', linewidths=0.4)

step_lbl = ax_f.text(0.97, 0.97, "", transform=ax_f.transAxes,
                     ha='right', va='top', fontsize=11,
                     fontfamily='monospace', color='white', alpha=0.6)
pop_lbl = ax_f.text(0.03, 0.97, "", transform=ax_f.transAxes,
                    ha='left', va='top', fontsize=10,
                    color='#aaccff', alpha=0.7)

# Stats
ax_s.set_facecolor('#080c18')
for s in ax_s.spines.values():
    s.set_color('#222244')
ax_s.tick_params(colors='#555555', labelsize=8)

steps_arr = [f['step'] for f in frames]
meanK_arr = [f['mean_K'] for f in frames]
n_arr = [f['n'] for f in frames]

ax_s.set_xlim(0, steps_arr[-1])
ax_s.set_ylim(0, 90)
ax_s.set_xlabel("Time step", color='#999999', fontsize=10)
ax_s.set_ylabel("Mean K", color='#ff9944', fontsize=10)
line_K, = ax_s.plot([], [], color='#ff9944', lw=2, label='Mean K')

ax_n = ax_s.twinx()
ax_n.set_ylim(0, 90)
ax_n.set_ylabel("Population", color='#44aaff', fontsize=10)
ax_n.tick_params(colors='#555555', labelsize=8)
ax_n.spines['right'].set_color('#222244')
line_n, = ax_n.plot([], [], color='#44aaff', lw=1.5, alpha=0.6, label='Pop')

ax_s.legend(handles=[line_K, line_n], loc='upper left', fontsize=9,
            facecolor='#080c18', edgecolor='#333355', labelcolor='#cccccc')

fps = 30


def update(idx):
    f = frames[idx]
    if len(f['xs']) == 0:
        scatter.set_offsets(np.empty((0, 2)))
        scatter.set_sizes([])
        scatter.set_array(np.array([]))
    else:
        scatter.set_offsets(np.column_stack([f['xs'], f['ys']]))
        scatter.set_sizes(r_to_s(f['rs']))
        scatter.set_array(np.clip((f['Ks'] - K_lo) / (K_hi - K_lo), 0, 1))

    step_lbl.set_text(f"t = {f['step']}")
    pop_lbl.set_text(f"n = {f['n']}   ⟨K⟩ = {f['mean_K']:.1f}")
    line_K.set_data(steps_arr[:idx + 1], meanK_arr[:idx + 1])
    line_n.set_data(steps_arr[:idx + 1], n_arr[:idx + 1])
    return [scatter, step_lbl, pop_lbl, line_K, line_n]


print("  Rendering MP4...")
ani = animation.FuncAnimation(fig, update, frames=len(frames),
                              interval=1000 / fps, blit=True)
writer = animation.FFMpegWriter(
    fps=fps, bitrate=3000,
    extra_args=['-vf', 'pad=ceil(iw/2)*2:ceil(ih/2)*2',
                '-pix_fmt', 'yuv420p'])

ani.save(str(OUT_MP4), writer=writer, dpi=120)
print(f"  saved {OUT_MP4}")
ani.save(str(OUT_REPO), writer=writer, dpi=120)
print(f"  saved {OUT_REPO}")
plt.close(fig)

import subprocess
for t_sec, label in [(0, "proto_start"), (6, "proto_mid"), (12, "proto_end")]:
    subprocess.run(["ffmpeg", "-y", "-ss", str(t_sec), "-i", str(OUT_REPO),
                    "-frames:v", "1", f"/tmp/verify_{label}.png"],
                   capture_output=True)
print("  Verification: /tmp/verify_proto_*.png")
print("Done.")
