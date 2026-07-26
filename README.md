# MOF-5 CO₂ Adsorption and Diffusion Simulation

This repository contains the complete computational workflow for studying CO₂ adsorption and self-diffusion in the metal–organic framework (MOF-5) using **RASPA2** and **LAMMPS**. The project combines Grand Canonical Monte Carlo (GCMC) simulations with Molecular Dynamics (MD) simulations to investigate CO₂ uptake and transport within the porous framework.

---

## Project Workflow

```
QMOF Database
      │
      ▼
Obtain MOF-5 CIF Structure
      │
      ▼
Generate Force Field Parameters
(lammps-interface + UFF4MOF)
      │
      ▼
RASPA2 GCMC Simulation
(CO₂ Adsorption)
      │
      ▼
Convert RASPA Output to LAMMPS
      │
      ▼
Energy Minimization
      │
      ▼
NVT Equilibration
      │
      ▼
NVE Production Run
      │
      ▼
CO₂ Diffusion Analysis
(MSD & Self-Diffusion Coefficient)
```

---

## Repository Structure

```
.
├── RAPSA2/                  # RASPA2 input and generated files
├── diffusion_analysis/      # Python scripts and MSD results
├── *.cif                    # MOF structure files
├── *.data                   # LAMMPS data files
├── *.in                     # LAMMPS input scripts
├── *.py                     # Analysis scripts
├── README.md
└── .gitignore
```

---

## Software Used

- RASPA2
- LAMMPS
- lammps-interface
- Python 3
- NumPy
- Matplotlib
- OVITO (trajectory visualization)

---

## Simulation Workflow

### 1. Framework Preparation

- Obtain the MOF-5 crystal structure from the QMOF database.
- Generate LAMMPS-compatible files using **lammps-interface** with the **UFF4MOF** force field.

### 2. GCMC Simulation (RASPA2)

- Temperature: **298 K**
- Pressure: **1 bar**
- Adsorbate: **CO₂ (TraPPE model)**

The GCMC simulation determines the equilibrium loading of CO₂ inside MOF-5.

### 3. Molecular Dynamics (LAMMPS)

The loaded structure is imported into LAMMPS for diffusion calculations.

Simulation sequence:

1. Energy Minimization
2. NVT Equilibration
3. NVE Production Simulation

The MOF framework is treated as rigid while CO₂ molecules evolve dynamically.

### 4. Post-processing

Mean Squared Displacement (MSD) is calculated using Python scripts to determine the self-diffusion coefficient of CO₂.

---

## Analysis

The repository includes scripts for:

- Reading LAMMPS trajectory files
- Computing Mean Squared Displacement (MSD)
- Linear fitting of the diffusive regime
- Calculating the CO₂ self-diffusion coefficient
- Plotting diffusion curves

---

## Notes

Large generated files such as:

- trajectory dumps (`*.dump`)
- restart files (`*.restart`)
- RASPA generated Movies/
- VTK files

are intentionally excluded from version control using `.gitignore`.

---

## Citation

If you use this repository in academic work, please cite the appropriate software:

- RASPA2
- LAMMPS
- UFF4MOF
- TraPPE CO₂ Force Field

---

## Author

**Md. Rakib Hasan**

Department of Electrical and Electronic Engineering (EEE)

Khulna University of Engineering & Technology (KUET)

Bangladesh

---

## License

This project is intended for academic and research purposes.
