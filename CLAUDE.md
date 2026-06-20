# Repo notes for Claude

## Purpose

Produce a Wolfram Community post — `community/biomolecular_condensates.nb` —
on the physics of biomolecular condensates (liquid-liquid phase separation),
their role as protocellular compartments, and connection to the origin of life.
First-principles mathematical modelling in both **Python** and **Wolfram Language**.

## Pipeline

### Wolfram Language (primary deliverable — the Community notebook)

```
wolfram/phase_diagram.wls       ── Flory-Huggins thermodynamics, spinodal/binodal
wolfram/cahn_hilliard.wls       ── 2D Cahn-Hilliard PDE simulation
wolfram/kinetics.wls            ── concentration enhancement, reaction kinetics
wolfram/protocell.wls           ── RNA replication, growth, division, selection
wolfram/disease.wls             ── liquid→gel→aggregate transition
wolfram/run_all.wls             ── one entry point, writes docs/images/*.png

community/build_notebook.wls    ── assembles biomolecular_condensates.nb + .pdf
```

### Python (parallel implementation in src/)

```
src/flory_huggins.py            ── thermodynamics, phase diagrams
src/cahn_hilliard.py            ── 2D PDE simulation with animation
src/condensate_kinetics.py      ── reaction rate enhancement model
src/protocell.py                ── RNA replication, growth, division, selection
src/requirements.txt            ── dependencies
```

Python figures go to `figures/`. Wolfram figures go to `docs/images/`.

## Conventions

* Plain-text `.wls` is the source of truth; the `.nb` and `.pdf` in
  `community/` are committed *outputs*.
* Figures live in `docs/images/` (Wolfram) and `figures/` (Python).
* All scripts are self-contained and runnable from the repo root.
* Wolfram execution: always via `/usr/local/bin/wolframscript`. Never
  use Wolfram MCP servers.

## Commit cadence

Commit + push after each meaningful step. Keep messages short and factual.
