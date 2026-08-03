# CO₂ Adsorption and Self-Diffusion in MOF-5 using GCMC–MD

A computational workflow for investigating CO₂ adsorption and self-diffusion in **MOF-5 (IRMOF-1)** using **RASPA2** and **LAMMPS**. The adsorption equilibrium is obtained from Grand Canonical Monte Carlo (GCMC) simulations, followed by Molecular Dynamics (MD) to evaluate the CO₂ self-diffusion coefficient.

---

## Overview

This project implements a multiscale **GCMC → MD** workflow:

1. Construct rigid MOF-5 framework.
2. Perform GCMC adsorption simulations at equilibrium.
3. Export the equilibrium configuration.
4. Equilibrate using NVT MD.
5. Perform NVE production MD.
6. Compute the CO₂ center-of-mass mean squared displacement (MSD).
7. Determine the self-diffusion coefficient using the Einstein relation.

---

## Software

- RASPA2
- LAMMPS (22 Jul 2025)
- Python 3
- NumPy
- Matplotlib

---

## System Description

| Property | Value |
|----------|-------|
| Material | MOF-5 (IRMOF-1) |
| Framework | Fully rigid |
| Framework force field | UFF |
| Adsorbate | CO₂ |
| CO₂ model | TraPPE (rigid) |
| Electrostatics | Ewald summation |
| Simulation cell | 3 × 3 × 2 supercell |
| Cell dimensions | 55.26 × 55.26 × 36.84 Å³ |

---

# Workflow

## Stage 1 — GCMC Adsorption (RASPA2)

Simulation conditions

- Ensemble: μVT (Grand Canonical Monte Carlo)
- Temperature: 298 K
- Pressure: 1 bar
- Initialization cycles: 10,000
- Production cycles: 50,000

### Results

| Property | Value |
|----------|------:|
| Average CO₂ loading | 19.27 molecules/unit cell |
| Final configuration | 19 molecules |
| Adsorption capacity | 30.59 mg g⁻¹ |
| Loading | 0.6952 mol kg⁻¹ |
| Host–CO₂ interaction energy | −24,717.20 K |
| CO₂–CO₂ interaction energy | −471.94 K |

The adsorption loading converged successfully, and the final equilibrium configuration containing **19 CO₂ molecules** was exported for Molecular Dynamics simulations.

---

## Stage 2 — NVT Equilibration (LAMMPS)

Simulation conditions

- Ensemble: NVT
- Framework: Fully frozen
- CO₂: TraPPE rigid molecules
- Time step: 1 fs

Purpose

- Relax the adsorbed configuration
- Thermalize the adsorbed CO₂ molecules before production dynamics

---

## Stage 3 — Production Diffusion (LAMMPS)

Simulation conditions

- Ensemble: NVE
- Framework: Frozen
- CO₂: Rigid (`fix rigid/nve/small`)
- Time step: 1 fs
- Production length: 5 ns

### Diffusion Results

| Property | Value |
|----------|------:|
| CO₂ molecules | 19 |
| Diffusive fit window | 500–1461 ps |
| Self-diffusion coefficient | **1.616 × 10⁻⁴ cm²/s** |
| Block averaged D | **(1.700 ± 0.178) × 10⁻⁴ cm²/s** |
| Dx | 1.825 × 10⁻⁴ cm²/s |
| Dy | 1.911 × 10⁻⁴ cm²/s |
| Dz | 1.112 × 10⁻⁴ cm²/s |
| Directional anisotropy | 49.4% |

The diffusion coefficient was calculated from the center-of-mass MSD using the Einstein relation.

---

## Analysis

The analysis script performs:

- Center-of-mass trajectory reconstruction
- MSD calculation
- Automatic identification of the Fickian diffusion regime
- Linear regression of the MSD
- Block averaging
- Per-axis diffusion analysis
- Per-molecule MSD analysis

Generated figures:

- `diffusive_regime_check.png`
- `msd_vs_time.png`
- `msd_per_axis.png`
- `msd_per_molecule.png`

---

## Repository Structure

```text
.
├── raspa/
│   ├── simulation.input
│   └── output/
│
├── lammps/
│   ├── in.nvt
│   ├── in.nve
│   └── restart/
│
├── analysis/
│   ├── diffusion_analysis.py
│   ├── msd_vs_time.png
│   ├── msd_per_axis.png
│   ├── msd_per_molecule.png
│   └── summary.txt
│
└── README.md
```

---

## Notes

- The MOF-5 framework was treated as completely rigid.
- Diffusion coefficients correspond to the equilibrium CO₂ loading obtained from GCMC.
- The reported diffusion coefficient is based on a single 5 ns trajectory.
- Longer production simulations or multiple independent trajectories are recommended for improved statistical convergence.

---

## Citation

If you use this workflow or repository, please cite the associated publication (to be added).

---

## Author

**Rakib Hasan**

Department of Electrical and Electronic Engineering (EEE)

Khulna University of Engineering & Technology (KUET)

Bangladesh
