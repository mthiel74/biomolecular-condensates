#!/usr/bin/env python3
"""Protocell dynamics: growth, division, competition, and selection.

Models the Lifshitz-Slyozov coarsening of condensate droplets and
their behaviour as proto-Darwinian entities. This is the full
origin-of-life argument: condensate droplets that concentrate
reactants, catalyse replication faster than degradation, grow by
accretion, divide by mechanical instability, and compete for
resources — primitive natural selection without genetic material.

Also models the liquid-to-gel-to-aggregate transition relevant to
neurodegenerative disease (FUS, TDP-43, tau protein condensates).
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from pathlib import Path

FIGURES = Path(__file__).resolve().parent.parent / "figures"
FIGURES.mkdir(exist_ok=True)


# =====================================================================
# 1. Lifshitz-Slyozov coarsening (Ostwald ripening)
# =====================================================================

def lifshitz_slyozov_radii(n_droplets=200, t_max=200.0, dt=0.01,
                           gamma=1.0, c_inf=1.0, D=1.0, Vm=1.0,
                           seed=42):
    """Simulate Ostwald ripening of a polydisperse droplet population.

    The Gibbs-Thomson equation gives the equilibrium concentration at
    a droplet of radius R:
        c_eq(R) = c_inf · exp(2γVm / RkT) ≈ c_inf (1 + l_c/R)

    where l_c = 2γVm/kT is the capillary length. Droplets with
    R < R_critical dissolve; droplets with R > R_critical grow.

    dR/dt = D Vm (c_∞ - c_eq(R)) / R = D Vm c_inf l_c (1/R* - 1/R) / R

    where R* is the critical radius (mean-field: R* = <R>).
    """
    rng = np.random.default_rng(seed)
    l_c = 2 * gamma * Vm  # capillary length (setting kT=1)

    # Initial radii: log-normal distribution
    R = rng.lognormal(mean=np.log(1.0), sigma=0.5, size=n_droplets)
    R = np.maximum(R, 0.01)

    times = [0.0]
    R_history = [R.copy()]
    n_alive = [n_droplets]

    t = 0.0
    while t < t_max:
        alive = R > 0.01
        if np.sum(alive) < 2:
            break

        R_star = np.mean(R[alive])

        dRdt = np.zeros_like(R)
        dRdt[alive] = D * Vm * c_inf * l_c * (1.0 / R_star - 1.0 / R[alive]) / R[alive]

        R += dRdt * dt
        R = np.maximum(R, 0.0)
        R[R < 0.01] = 0.0  # dissolve tiny droplets

        t += dt
        if int(t / dt) % max(1, int(1.0 / dt)) == 0:
            times.append(t)
            R_history.append(R.copy())
            n_alive.append(np.sum(R > 0.01))

    return np.array(times), R_history, np.array(n_alive)


# =====================================================================
# 2. Protocell population dynamics with division and selection
# =====================================================================

class Protocell:
    """A single protocell (condensate droplet) with internal chemistry."""

    def __init__(self, radius, K, replication_rate, rng):
        self.radius = radius
        self.K = K
        self.replication_rate = replication_rate
        self.volume = (4/3) * np.pi * radius**3
        self.internal_concentration = 0.0
        self.age = 0.0
        self.generation = 0
        self.rng = rng

    def grow(self, c_bulk, dt):
        """Grow by accretion of material from the bulk."""
        self.internal_concentration = self.K * c_bulk
        growth_rate = self.replication_rate * self.internal_concentration
        degradation_rate = 0.05 * self.internal_concentration

        net_rate = growth_rate - degradation_rate
        if net_rate > 0:
            dV = net_rate * self.volume * dt * 0.01
            self.volume += dV
            self.radius = (3 * self.volume / (4 * np.pi))**(1/3)

        self.age += dt

    def should_divide(self, R_div=3.0):
        """Division when radius exceeds threshold."""
        return self.radius > R_div

    def divide(self):
        """Split into two daughter cells with slight variation."""
        r_new = self.radius / 2**(1/3)
        K_noise = 1.0 + 0.05 * self.rng.standard_normal()
        k_noise = 1.0 + 0.05 * self.rng.standard_normal()

        d1 = Protocell(r_new, self.K * K_noise,
                       self.replication_rate * k_noise, self.rng)
        d2 = Protocell(r_new, self.K / K_noise,
                       self.replication_rate / k_noise, self.rng)
        d1.generation = self.generation + 1
        d2.generation = self.generation + 1
        return d1, d2


def run_protocell_simulation(n_initial=50, t_max=500.0, dt=0.5,
                             c_bulk=0.1, seed=42):
    """Simulate a population of competing protocells."""
    rng = np.random.default_rng(seed)

    cells = []
    for _ in range(n_initial):
        r = rng.uniform(0.5, 2.0)
        K = rng.uniform(10, 200)
        k_rep = rng.uniform(0.01, 0.2)
        cells.append(Protocell(r, K, k_rep, rng))

    times = []
    n_cells = []
    mean_K = []
    mean_k_rep = []
    mean_radius = []
    max_gen = []

    t = 0.0
    while t < t_max and len(cells) > 0:
        # Grow
        for cell in cells:
            cell.grow(c_bulk, dt)

        # Division
        new_cells = []
        for cell in cells:
            if cell.should_divide():
                d1, d2 = cell.divide()
                new_cells.extend([d1, d2])
            else:
                new_cells.append(cell)
        cells = new_cells

        # Death: remove cells with R < 0.1
        cells = [c for c in cells if c.radius > 0.1]

        # Resource competition: limit population
        max_pop = 500
        if len(cells) > max_pop:
            fitness = np.array([c.K * c.replication_rate for c in cells])
            keep_idx = np.argsort(fitness)[-max_pop:]
            cells = [cells[i] for i in keep_idx]

        # Record
        if int(t / dt) % 10 == 0 and len(cells) > 0:
            times.append(t)
            n_cells.append(len(cells))
            mean_K.append(np.mean([c.K for c in cells]))
            mean_k_rep.append(np.mean([c.replication_rate for c in cells]))
            mean_radius.append(np.mean([c.radius for c in cells]))
            max_gen.append(max(c.generation for c in cells))

        t += dt

    return {
        'times': np.array(times),
        'n_cells': np.array(n_cells),
        'mean_K': np.array(mean_K),
        'mean_k_rep': np.array(mean_k_rep),
        'mean_radius': np.array(mean_radius),
        'max_gen': np.array(max_gen),
        'final_cells': cells
    }


# =====================================================================
# 3. Liquid → gel → aggregate transition (disease)
# =====================================================================

def aging_transition_model(t_max=100.0, dt=0.1, k_gel=0.02, k_agg=0.005):
    """Model the liquid → gel → solid aggregate transition.

    Three-state kinetic model:
        d[liquid]/dt = -k_gel · [liquid]
        d[gel]/dt    = k_gel · [liquid] - k_agg · [gel]
        d[aggregate]/dt = k_agg · [gel]
    """
    n_steps = int(t_max / dt)
    t = np.zeros(n_steps)
    liquid = np.zeros(n_steps)
    gel = np.zeros(n_steps)
    aggregate = np.zeros(n_steps)

    liquid[0] = 1.0

    for i in range(1, n_steps):
        t[i] = t[i-1] + dt
        liquid[i] = liquid[i-1] - k_gel * liquid[i-1] * dt
        gel[i] = gel[i-1] + (k_gel * liquid[i-1] - k_agg * gel[i-1]) * dt
        aggregate[i] = aggregate[i-1] + k_agg * gel[i-1] * dt

    return t, liquid, gel, aggregate


# =====================================================================
# Plotting
# =====================================================================

def plot_ostwald_ripening():
    """Visualise Ostwald ripening of droplet population."""
    times, R_history, n_alive = lifshitz_slyozov_radii(
        n_droplets=150, t_max=150.0, dt=0.01)

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    # Panel 1: size distribution at several times
    ax = axes[0, 0]
    time_indices = [0, len(times)//4, len(times)//2, -1]
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(time_indices)))
    for idx, c in zip(time_indices, colors):
        R = R_history[idx]
        R_alive = R[R > 0.01]
        if len(R_alive) > 0:
            ax.hist(R_alive, bins=20, alpha=0.5, color=c,
                    label=f't = {times[idx]:.0f}', density=True)
    ax.set_xlabel("Droplet radius R", fontsize=12)
    ax.set_ylabel("Probability density", fontsize=12)
    ax.set_title("Size distribution evolution", fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    # Panel 2: number of droplets vs time
    ax = axes[0, 1]
    ax.plot(times, n_alive, 'b-', lw=2)
    ax.set_xlabel("Time", fontsize=12)
    ax.set_ylabel("Number of surviving droplets", fontsize=12)
    ax.set_title("Droplet count: Ostwald ripening", fontsize=13)
    ax.grid(True, alpha=0.3)

    # Panel 3: mean and max radius vs time
    ax = axes[1, 0]
    mean_R = [np.mean(R[R > 0.01]) if np.any(R > 0.01) else 0 for R in R_history]
    max_R = [np.max(R) if np.any(R > 0.01) else 0 for R in R_history]
    ax.plot(times, mean_R, 'b-', lw=2, label='Mean radius')
    ax.plot(times, max_R, 'r--', lw=2, label='Max radius')
    t_fit = times[len(times)//3:]
    if len(t_fit) > 0 and t_fit[0] > 0:
        c_fit = mean_R[len(times)//3] / times[len(times)//3]**(1/3)
        ax.plot(t_fit, c_fit * t_fit**(1/3), 'k:', lw=1.5,
                label=r'$\langle R \rangle \sim t^{1/3}$ (LSW theory)')
    ax.set_xlabel("Time", fontsize=12)
    ax.set_ylabel("Radius", fontsize=12)
    ax.set_title("Coarsening: Lifshitz–Slyozov–Wagner", fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    # Panel 4: schematic of the process
    ax = axes[1, 1]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_aspect('equal')
    ax.set_title("Ostwald ripening schematic", fontsize=13)

    # Small droplets (dissolving)
    for x, y, r in [(2, 3, 0.3), (3, 7, 0.2), (1.5, 6, 0.15),
                     (4, 2, 0.25), (2.5, 5, 0.1)]:
        circle = Circle((x, y), r, fill=True, fc='lightsalmon',
                        ec='red', lw=1, alpha=0.6)
        ax.add_patch(circle)
        ax.annotate('', xy=(5, 5), xytext=(x, y),
                    arrowprops=dict(arrowstyle='->', color='red', alpha=0.3, lw=0.8))

    # Large droplet (growing)
    big = Circle((6.5, 5), 1.8, fill=True, fc='lightblue',
                 ec='navy', lw=2.5, alpha=0.7)
    ax.add_patch(big)
    ax.text(6.5, 5, 'Growing\ndroplet', ha='center', va='center',
            fontsize=10, fontweight='bold', color='navy')

    ax.text(2, 1, 'Small droplets\ndissolve', ha='center',
            fontsize=10, color='red', fontstyle='italic')
    ax.text(6.5, 1, 'Large droplets\ngrow', ha='center',
            fontsize=10, color='navy', fontstyle='italic')
    ax.axis('off')

    fig.suptitle("Lifshitz–Slyozov coarsening of condensate droplets",
                 fontsize=15, y=1.01)
    fig.tight_layout()
    fig.savefig(FIGURES / "ostwald_ripening.png", dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"  saved {FIGURES / 'ostwald_ripening.png'}")


def plot_protocell_selection():
    """Show proto-Darwinian selection on condensate populations."""
    print("  Running protocell simulation...")
    result = run_protocell_simulation(n_initial=50, t_max=500.0)

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    # Panel 1: population size
    ax = axes[0, 0]
    ax.plot(result['times'], result['n_cells'], 'b-', lw=2)
    ax.set_xlabel("Time", fontsize=12)
    ax.set_ylabel("Population size", fontsize=12)
    ax.set_title("Protocell population dynamics", fontsize=13)
    ax.grid(True, alpha=0.3)

    # Panel 2: mean partition coefficient (should increase via selection)
    ax = axes[0, 1]
    ax.plot(result['times'], result['mean_K'], 'r-', lw=2)
    ax.set_xlabel("Time", fontsize=12)
    ax.set_ylabel("Mean partition coefficient K", fontsize=12)
    ax.set_title("Selection for higher K", fontsize=13)
    ax.grid(True, alpha=0.3)

    # Panel 3: mean replication rate (should increase)
    ax = axes[1, 0]
    ax.plot(result['times'], result['mean_k_rep'], 'g-', lw=2)
    ax.set_xlabel("Time", fontsize=12)
    ax.set_ylabel("Mean replication rate", fontsize=12)
    ax.set_title("Selection for faster replication", fontsize=13)
    ax.grid(True, alpha=0.3)

    # Panel 4: generation count
    ax = axes[1, 1]
    ax.plot(result['times'], result['max_gen'], 'purple', lw=2)
    ax.set_xlabel("Time", fontsize=12)
    ax.set_ylabel("Max generation", fontsize=12)
    ax.set_title("Lineage depth", fontsize=13)
    ax.grid(True, alpha=0.3)

    fig.suptitle("Proto-Darwinian selection in condensate populations",
                 fontsize=15, y=1.01)
    fig.tight_layout()
    fig.savefig(FIGURES / "protocell_selection.png", dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"  saved {FIGURES / 'protocell_selection.png'}")


def plot_disease_transition():
    """Plot the liquid → gel → aggregate transition."""
    t1, liq1, gel1, agg1 = aging_transition_model(k_gel=0.02, k_agg=0.005)
    t2, liq2, gel2, agg2 = aging_transition_model(k_gel=0.05, k_agg=0.01)
    t3, liq3, gel3, agg3 = aging_transition_model(k_gel=0.1, k_agg=0.02)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    for ax, (t, liq, gel, agg, title) in zip(axes, [
        (t1, liq1, gel1, agg1, "Slow maturation\n(k_gel=0.02, k_agg=0.005)"),
        (t2, liq2, gel2, agg2, "Moderate\n(k_gel=0.05, k_agg=0.01)"),
        (t3, liq3, gel3, agg3, "Fast (disease mutant)\n(k_gel=0.1, k_agg=0.02)")
    ]):
        ax.fill_between(t, 0, liq, alpha=0.3, color='dodgerblue', label='Liquid')
        ax.fill_between(t, liq, liq + gel, alpha=0.3, color='orange', label='Gel')
        ax.fill_between(t, liq + gel, liq + gel + agg, alpha=0.3, color='red', label='Aggregate')
        ax.plot(t, liq, 'b-', lw=2)
        ax.plot(t, gel, color='orange', lw=2)
        ax.plot(t, agg, 'r-', lw=2)
        ax.set_xlabel("Time (a.u.)", fontsize=12)
        ax.set_ylabel("Fraction", fontsize=12)
        ax.set_title(title, fontsize=12)
        ax.legend(fontsize=9, loc='right')
        ax.set_ylim(0, 1.05)
        ax.grid(True, alpha=0.3)

    fig.suptitle("Condensate aging: liquid → gel → aggregate\n"
                 "(FUS, TDP-43, tau in neurodegeneration)",
                 fontsize=14, y=1.04)
    fig.tight_layout()
    fig.savefig(FIGURES / "disease_transition.png", dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"  saved {FIGURES / 'disease_transition.png'}")


def plot_disease_phase_diagram():
    """Phase diagram of condensate material states."""
    fig, ax = plt.subplots(figsize=(8, 6))

    # Axes: temperature (or molecular interactions) vs concentration
    conc = np.linspace(0, 10, 200)
    T = np.linspace(0, 10, 200)
    CC, TT = np.meshgrid(conc, T)

    # Schematic boundaries
    # Liquid-liquid phase boundary (binodal-like)
    boundary_ll = 3.0 + 2.0 * np.exp(-0.5 * conc)
    # Gel boundary
    boundary_gel = 2.0 + 1.0 * np.exp(-0.3 * conc)
    # Aggregate boundary
    boundary_agg = 1.0 + 0.5 * np.exp(-0.2 * conc)

    ax.fill_between(conc, boundary_ll, 10, alpha=0.15, color='dodgerblue')
    ax.fill_between(conc, boundary_gel, boundary_ll, alpha=0.15, color='orange')
    ax.fill_between(conc, boundary_agg, boundary_gel, alpha=0.15, color='red')
    ax.fill_between(conc, 0, boundary_agg, alpha=0.15, color='darkred')

    ax.plot(conc, boundary_ll, 'b-', lw=2, label='Liquid condensate boundary')
    ax.plot(conc, boundary_gel, color='orange', lw=2, label='Gel transition')
    ax.plot(conc, boundary_agg, 'r-', lw=2, label='Solid aggregate')

    ax.text(5, 8, 'Mixed\n(one phase)', fontsize=13, ha='center',
            color='dodgerblue', fontweight='bold')
    ax.text(7, 4.5, 'Liquid\ncondensate', fontsize=13, ha='center',
            color='darkorange', fontweight='bold')
    ax.text(7, 2.5, 'Gel', fontsize=12, ha='center',
            color='orangered', fontweight='bold')
    ax.text(7, 0.7, 'Solid\naggregate', fontsize=12, ha='center',
            color='darkred', fontweight='bold')

    ax.annotate('Disease mutations\nshift boundaries →',
                xy=(3, 3), fontsize=10, color='red',
                xytext=(0.5, 4.5),
                arrowprops=dict(arrowstyle='->', color='red', lw=2))

    ax.set_xlabel("Protein concentration (a.u.)", fontsize=13)
    ax.set_ylabel("Temperature / Interaction strength (a.u.)", fontsize=13)
    ax.set_title("Condensate material state phase diagram", fontsize=14)
    ax.legend(fontsize=10, loc='upper right')
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(FIGURES / "disease_phase_diagram.png", dpi=200)
    plt.close(fig)
    print(f"  saved {FIGURES / 'disease_phase_diagram.png'}")


if __name__ == "__main__":
    print("Protocell dynamics and disease transitions")
    print("=" * 40)

    print("\nGenerating figures...")
    plot_ostwald_ripening()
    plot_protocell_selection()
    plot_disease_transition()
    plot_disease_phase_diagram()
    print("Done.")
