# Biomolecular Condensates and the Origin of Life

First-principles mathematical modelling of liquid-liquid phase separation
in biology — from Flory-Huggins thermodynamics to protocellular compartments.

## What this is

Biomolecular condensates are membrane-free liquid droplets inside cells, formed
by liquid-liquid phase separation (LLPS) of proteins and nucleic acids.
Discovered as P granules in *C. elegans* germ cells (Brangwynne et al., 2009),
they organise gene regulation, stress response, and signal transduction. They
also revive Oparin's 1924 coacervate hypothesis for the origin of life: phase
separation concentrates dilute reactants by factors of 100-1000, solving the
concentration problem in prebiotic chemistry.

This repository implements the mathematical physics of condensates in both
**Python** and **Wolfram Language**, covering six topics:

1. **Flory-Huggins thermodynamics** — free energy, spinodal, binodal, critical point
2. **Cahn-Hilliard dynamics** — spinodal decomposition PDE (2D spectral in Python, 1D Method of Lines in WL)
3. **Concentration enhancement** — partition coefficient K, bimolecular rate enhancement K^2
4. **RNA world kinetics** — Michaelis-Menten ribozyme replication vs hydrolytic degradation
5. **Protocell dynamics** — Ostwald ripening, growth/division, proto-Darwinian selection
6. **Disease transitions** — liquid to gel to aggregate (FUS, TDP-43, tau in neurodegeneration)

## Repository structure

```
src/
  flory_huggins.py          Thermodynamics, phase diagrams
  cahn_hilliard.py          2D semi-implicit spectral solver
  condensate_kinetics.py    Reaction rate enhancement model
  protocell.py              Ostwald ripening, protocells, disease model

wolfram/
  phase_diagram.wls         Flory-Huggins in WL
  cahn_hilliard.wls         1D Cahn-Hilliard via NDSolve
  kinetics.wls              Rate enhancement, ribozyme kinetics
  protocell.wls             Ostwald ripening, protocell selection
  disease.wls               Liquid-gel-aggregate transition
  run_all.wls               Master runner for all WL scripts

community/
  build_notebook.wls        Notebook builder (ENSO-emergence pattern)
  biomolecular_condensates.nb   Generated Wolfram Community notebook
  biomolecular_condensates.pdf  PDF export

figures/                    Python-generated figures (15 PNGs)
docs/images/                Wolfram-generated figures (13 PNGs)
```

## Quick start

### Python

```bash
pip install numpy scipy matplotlib
python3 src/flory_huggins.py
python3 src/cahn_hilliard.py
python3 src/condensate_kinetics.py
python3 src/protocell.py
```

### Wolfram Language

```bash
wolframscript -file wolfram/run_all.wls        # all analysis scripts
wolframscript -file community/build_notebook.wls  # rebuild the notebook
```

## Key equations

**Flory-Huggins free energy** (per lattice site):

$$\frac{f}{kT} = \frac{\phi \ln \phi}{N_1} + \frac{(1-\phi)\ln(1-\phi)}{N_2} + \chi\,\phi\,(1-\phi)$$

**Cahn-Hilliard equation** (conserved dynamics):

$$\frac{\partial \phi}{\partial t} = M \nabla^2 \left(\frac{\delta F}{\delta \phi} - \kappa \nabla^2 \phi \right)$$

**Bimolecular rate enhancement**:

$$\text{rate}_{\text{condensate}} = K_A \, K_B \, \text{rate}_{\text{dilute}}$$

## Key references

- Brangwynne, C. P. et al. (2009). *Science* 324, 1729-1732.
- Hyman, A. A., Weber, C. A., & Julicher, F. (2014). *Annu. Rev. Cell Dev. Biol.* 30, 39-58.
- Banani, S. F. et al. (2017). *Nat. Rev. Mol. Cell Biol.* 18, 285-298.
- Shin, Y. & Brangwynne, C. P. (2017). *Science* 357, eaaf4382.
- Patel, A. et al. (2015). *Cell* 162, 1066-1077.
- Zwicker, D. et al. (2017). *Nat. Phys.* 13, 408-413.

## License

MIT
