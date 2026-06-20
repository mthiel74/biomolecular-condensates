#!/usr/bin/env python3
"""2D Cahn-Hilliard simulation of spinodal decomposition.

The Cahn-Hilliard equation governs the time evolution of a conserved
order parameter φ(x,t) during phase separation:

    ∂φ/∂t = M ∇²[ ∂f/∂φ − κ ∇²φ ]

where f(φ) is the Flory-Huggins free energy density, M is the mobility,
and κ is the gradient energy coefficient (surface tension parameter).

We solve this on a 2D periodic domain using a semi-implicit spectral
method: the biharmonic (linear) part is treated implicitly for stability,
while the nonlinear chemical potential is explicit.
"""

import numpy as np
from numpy.fft import fft2, ifft2, fftfreq
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from pathlib import Path

FIGURES = Path(__file__).resolve().parent.parent / "figures"
FIGURES.mkdir(exist_ok=True)


def flory_huggins_mu(phi, chi, N1=100, N2=1):
    """Chemical potential df/dφ from Flory-Huggins theory."""
    eps = 1e-10
    phi_safe = np.clip(phi, eps, 1 - eps)
    return ((np.log(phi_safe) + 1) / N1
            - (np.log(1 - phi_safe) + 1) / N2
            + chi * (1 - 2 * phi_safe))


class CahnHilliard2D:
    """Semi-implicit spectral solver for the 2D Cahn-Hilliard equation.

    The scheme in Fourier space:

        φ̂^{n+1} = (φ̂^n + dt M k² μ̂^n) / (1 + dt M κ k⁴)

    where k² = kx² + ky² is the squared wavenumber.
    """

    def __init__(self, L=10.0, N=128, chi=1.5, N1=100, N2=1,
                 mobility=1.0, kappa=0.5, dt=0.005):
        self.L = L
        self.N = N
        self.chi = chi
        self.N1 = N1
        self.N2 = N2
        self.M = mobility
        self.kappa = kappa
        self.dt = dt
        self.dx = L / N

        kx = 2 * np.pi * fftfreq(N, d=self.dx)
        ky = 2 * np.pi * fftfreq(N, d=self.dx)
        KX, KY = np.meshgrid(kx, ky)
        self.k2 = KX**2 + KY**2
        self.k4 = self.k2**2

        self.phi = None
        self.time = 0.0
        self.step_count = 0

    def initialize(self, phi_mean=0.3, noise_amplitude=0.01, seed=42):
        """Random initial condition near mean composition."""
        rng = np.random.default_rng(seed)
        self.phi = phi_mean + noise_amplitude * rng.standard_normal((self.N, self.N))
        self.phi = np.clip(self.phi, 1e-6, 1 - 1e-6)
        self.time = 0.0
        self.step_count = 0

    def step(self):
        """One semi-implicit time step."""
        mu = flory_huggins_mu(self.phi, self.chi, self.N1, self.N2)
        mu_hat = fft2(mu)
        phi_hat = fft2(self.phi)

        numerator = phi_hat + self.dt * self.M * self.k2 * mu_hat
        denominator = 1.0 + self.dt * self.M * self.kappa * self.k4
        phi_hat_new = numerator / denominator

        self.phi = np.real(ifft2(phi_hat_new))
        self.phi = np.clip(self.phi, 1e-6, 1 - 1e-6)
        self.time += self.dt
        self.step_count += 1

    def run(self, n_steps):
        """Run n_steps time steps."""
        for _ in range(n_steps):
            self.step()

    def total_free_energy(self):
        """Compute the total Flory-Huggins + gradient free energy."""
        eps = 1e-10
        phi = np.clip(self.phi, eps, 1 - eps)
        f_local = (phi * np.log(phi) / self.N1
                   + (1 - phi) * np.log(1 - phi) / self.N2
                   + self.chi * phi * (1 - phi))

        grad_x = np.roll(phi, -1, axis=1) - np.roll(phi, 1, axis=1)
        grad_y = np.roll(phi, -1, axis=0) - np.roll(phi, 1, axis=0)
        grad_sq = (grad_x / (2 * self.dx))**2 + (grad_y / (2 * self.dx))**2
        f_grad = 0.5 * self.kappa * grad_sq

        return np.sum(f_local + f_grad) * self.dx**2

    def droplet_count(self, threshold=0.5):
        """Rough count of dense-phase domains above threshold."""
        from scipy.ndimage import label
        binary = self.phi > threshold
        _, n_features = label(binary)
        return n_features

    def domain_size(self):
        """Characteristic domain size from structure factor peak."""
        phi_fluct = self.phi - np.mean(self.phi)
        S = np.abs(fft2(phi_fluct))**2
        k_mag = np.sqrt(self.k2)
        k_bins = np.linspace(0, np.max(k_mag) / 2, 50)
        S_radial = np.zeros(len(k_bins) - 1)
        for i in range(len(k_bins) - 1):
            mask = (k_mag >= k_bins[i]) & (k_mag < k_bins[i+1])
            if np.any(mask):
                S_radial[i] = np.mean(S[mask])
        k_centers = 0.5 * (k_bins[:-1] + k_bins[1:])
        if np.max(S_radial[1:]) > 0:
            k_peak = k_centers[1 + np.argmax(S_radial[1:])]
            return 2 * np.pi / k_peak if k_peak > 0 else self.L
        return self.L


def run_spinodal_decomposition(total_time=20.0, save_interval=None,
                               N=128, chi=1.5, phi_mean=0.3):
    """Run a full spinodal decomposition simulation and return snapshots."""
    sim = CahnHilliard2D(L=10.0, N=N, chi=chi, dt=0.005, kappa=0.5)
    sim.initialize(phi_mean=phi_mean)

    if save_interval is None:
        save_interval = total_time / 20

    n_total = int(total_time / sim.dt)
    n_save = max(1, int(save_interval / sim.dt))

    snapshots = [(0.0, sim.phi.copy(), sim.total_free_energy())]

    for i in range(1, n_total + 1):
        sim.step()
        if i % n_save == 0:
            snapshots.append((sim.time, sim.phi.copy(), sim.total_free_energy()))

    return sim, snapshots


def plot_spinodal_snapshots():
    """Generate multi-panel figure showing spinodal decomposition over time."""
    print("  Running Cahn-Hilliard simulation (this takes ~30s)...")
    sim, snapshots = run_spinodal_decomposition(
        total_time=50.0, save_interval=50.0 / 8, N=128, chi=1.5, phi_mean=0.3)

    times_to_show = [0.0, 2.0, 5.0, 10.0, 20.0, 50.0]
    panels = []
    for t_target in times_to_show:
        best = min(snapshots, key=lambda s: abs(s[0] - t_target))
        panels.append(best)

    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    axes = axes.ravel()

    for ax, (t, phi, F) in zip(axes, panels):
        im = ax.imshow(phi, extent=[0, sim.L, 0, sim.L],
                       cmap='RdBu_r', vmin=0, vmax=0.8,
                       origin='lower', interpolation='bilinear')
        ax.set_title(f"t = {t:.1f}", fontsize=13, fontweight='bold')
        ax.set_xlabel("x")
        ax.set_ylabel("y")

    fig.suptitle("Spinodal decomposition via Cahn–Hilliard equation\n"
                 f"(χ = {sim.chi}, φ̄ = 0.3, N = {sim.N1})",
                 fontsize=15, y=0.98)
    cbar = fig.colorbar(im, ax=axes, shrink=0.6, label='Volume fraction φ')
    fig.tight_layout(rect=[0, 0, 0.92, 0.94])
    fig.savefig(FIGURES / "cahn_hilliard_snapshots.png", dpi=200)
    plt.close(fig)
    print(f"  saved {FIGURES / 'cahn_hilliard_snapshots.png'}")

    return snapshots


def plot_free_energy_decay(snapshots):
    """Plot the total free energy vs time — should monotonically decrease."""
    times = [s[0] for s in snapshots]
    energies = [s[2] for s in snapshots]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(times, energies, 'b-', lw=2)
    ax.set_xlabel("Time", fontsize=13)
    ax.set_ylabel("Total free energy  F[φ]", fontsize=13)
    ax.set_title("Free energy decay during spinodal decomposition", fontsize=14)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGURES / "cahn_hilliard_energy.png", dpi=200)
    plt.close(fig)
    print(f"  saved {FIGURES / 'cahn_hilliard_energy.png'}")


def plot_coarsening(total_time=100.0, N=128):
    """Track domain size growth over time — should follow t^(1/3) Lifshitz-Slyozov."""
    print("  Running coarsening simulation (this takes ~60s)...")
    sim = CahnHilliard2D(L=20.0, N=N, chi=1.5, dt=0.005, kappa=0.5)
    sim.initialize(phi_mean=0.3)

    n_total = int(total_time / sim.dt)
    n_measure = max(1, n_total // 100)

    times = []
    domain_sizes = []

    for i in range(n_total):
        sim.step()
        if i % n_measure == 0 and sim.time > 1.0:
            times.append(sim.time)
            domain_sizes.append(sim.domain_size())

    times = np.array(times)
    domain_sizes = np.array(domain_sizes)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.loglog(times, domain_sizes, 'bo-', ms=3, lw=1.5, label='Simulation')

    t_fit = np.linspace(times[len(times)//3], times[-1], 100)
    c = domain_sizes[len(times)//2] / times[len(times)//2]**(1/3)
    ax.loglog(t_fit, c * t_fit**(1/3), 'r--', lw=2,
              label=r'$\ell \sim t^{1/3}$ (Lifshitz–Slyozov)')

    ax.set_xlabel("Time", fontsize=13)
    ax.set_ylabel("Domain size  ℓ", fontsize=13)
    ax.set_title("Coarsening dynamics: Ostwald ripening", fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, which='both')
    fig.tight_layout()
    fig.savefig(FIGURES / "cahn_hilliard_coarsening.png", dpi=200)
    plt.close(fig)
    print(f"  saved {FIGURES / 'cahn_hilliard_coarsening.png'}")


if __name__ == "__main__":
    print("Cahn-Hilliard spinodal decomposition")
    print("=" * 40)

    snapshots = plot_spinodal_snapshots()
    plot_free_energy_decay(snapshots)
    plot_coarsening(total_time=80.0)
    print("Done.")
