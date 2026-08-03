# Comprehensive Literature Comparison: CO₂ Adsorption & Self-Diffusion in MOF-5

## Executive Summary

Your simulation results are **physically reasonable and largely consistent with published Q1 literature** for CO₂ in MOF-5 (IRMOF-1). Your self-diffusion coefficient of **1.616 × 10⁻⁸ m²/s** (1.616 × 10⁻⁴ cm²/s) sits squarely within the narrow band of published MD self-diffusion coefficients (1.4–3.0 × 10⁻⁸ m²/s). Your adsorption capacity of **30.6 mg/g** at 298 K / 1 bar is also consistent with low-pressure GCMC and experimental data. However, there are **methodological red flags**—most notably the pronounced diffusion anisotropy in a cubic crystal—that require attention before journal submission.
I'll conduct a comprehensive literature search to find Q1 papers matching your criteria, then perform the detailed comparison and analysis you've requested.


---

# Comprehensive Literature Comparison: CO₂ Adsorption & Self-Diffusion in MOF-5

## Executive Summary

Your simulation results are **physically reasonable and largely consistent with published Q1 literature** for CO₂ in MOF-5 (IRMOF-1). Your self-diffusion coefficient of **1.616 × 10⁻⁸ m²/s** (1.616 × 10⁻⁴ cm²/s) sits squarely within the narrow band of published MD self-diffusion coefficients (1.4–3.0 × 10⁻⁸ m²/s). Your adsorption capacity of **30.6 mg/g** at 298 K / 1 bar is also consistent with low-pressure GCMC and experimental data. However, there are **methodological red flags**—most notably the pronounced diffusion anisotropy in a cubic crystal—that require attention before journal submission.

---

## Detailed Comparison Table

| Parameter | **Your Work** | Babarao & Jiang (2008) [1] | Listyarini *et al.* (2023) [2] | Zhao *et al.* (2009) [3] | Liu *et al.* (2025) [4] | Cheng *et al.* (2024) [5] | Saha *et al.* (2010) [6] | Choi *et al.* (2008) [7] |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| **Journal** | — | *Langmuir* | *J. Phys. Chem. B* | *Ind. Eng. Chem. Res.* | *Sep. Purif. Technol.* | *J. Membr. Sci.* | *Environ. Sci. Technol.* | *Microporous Mesoporous Mater.* |
| **Quartile** | — | Q1 | Q1 | Q1 | Q1 | Q1 | Q1 | Q1 |
| **Year** | 2026 | 2008 | 2023 | 2009 | 2025 | 2024 | 2010 | 2008 |
| **Force Field** | UFF + TraPPE | Classical FF (DREIDING/UFF-like) | SCC-DFTB/3ob/D3 | — (Experimental) | DFT + MD | UFF / TraPPE | — (Experimental) | — (Experimental) |
| **Framework** | Rigid | Rigid | Flexible (NPT) | Crystalline | Rigid/flexible | Rigid | Crystalline | Crystalline |
| **Temperature** | 298 K | 298–398 K | 248–398 K | 295–331 K | 298 K | 298 K | 298 K | 298 K |
| **Pressure** | 1 bar | Various | 1.013 bar | ≤1 atm | 0.1–0.3 MPa | Various | ≤1.05 bar | 40 bar |
| **CO₂ Model** | TraPPE rigid | — | DFTB-derived | — | — | TraPPE | — | — |
| **Supercell** | 3×3×2 (55.3×55.3×36.8 Å) | — | 1×1×1 (~26.1 Å) | 40–60 μm crystals | 3×3×2 (8.0×8.0×5.6 nm) | ≥24 Å cutoff | — | — |
| **Simulation Method** | GCMC → NVT → NVE MD | MD (NVT/NVE) | NPT MD | Gravimetric | GCMC + EMD + DFT | GCMC + MD | Volumetric/Gravimetric | Gravimetric |
| **Production Length** | 5 ns | — | 0.53 ns / temp | — | — | 5 ns (last) | — | — |
| **Loading (molecules)** | 19 in supercell | Various | 1, 4, 8, 12, 16 | — | Various | Various | — | — |
| **Loading (mol/kg)** | 0.695 | — | — | — | — | — | — | — |
| **Self-Diffusion Coefficient (MD)** | **1.616 × 10⁻⁴ cm²/s** (1.616 × 10⁻⁸ m²/s) | **1.4–3.0 × 10⁻⁸ m²/s** [cited in 2] | **2.09 × 10⁻⁸ m²/s** (1 CO₂) | — | — | Self-diffusivity reported | — | — |
| **Experimental Diffusion** | — | — | — | **8.1–11.5 × 10⁻⁹ cm²/s** (micropore) | — | — | **~10⁻⁹ m²/s** (uptake kinetics) | — |
| **Adsorption Capacity** | 30.59 mg/g (0.695 mmol/g) | — | — | — | — | Consistent with exp. | ~20–35 mg/g (est. at 1 bar) | ~812 mg/g (at 40 bar) |
| **Heat of Adsorption / Interaction Energy** | Host–CO₂: –10.8 kJ/mol (per molecule) | — | –13.0 kJ/mol (int. energy); exp. –15.1 kJ/mol | ~34 kJ/mol (isosteric) | — | Isosteric heat reported | — | 15.8–16.5 kJ/mol (low P) |
| **Activation Energy** | — | 4.05 kJ/mol | 5.47 kJ/mol | 7.61 kJ/mol | — | Reported | — | — |
| **MSD Analysis** | COM-MSD, Einstein, auto diffusive window | MSD, Einstein | MSD, Einstein | — | MSD | MSD | — | — |
| **Anisotropy Reported** | **Dx=1.825, Dy=1.911, Dz=1.112 (49.4%)** | Isotropic (cubic) | Isotropic (cubic) | — | — | — | — | — |
| **Agreement with Your Work** | — | **Excellent** (Ds within 15%) | **Very Good** (Ds within 23%) | N/A (different meas.) | **Methodology match** (same supercell) | Methodology match | Capacity consistent | Qst consistent |

---

### References

[1] R. Babarao, Z. Hu, J. Jiang, **"Diffusion and Separation of CO₂ and CH₄ in Silicalite, C168 Schwarzite, and IRMOF-1: A Comparative Study from Molecular Dynamics Simulation,"** *Langmuir*, 2008, 24(10), 5474–5484. DOI: [10.1021/la703434s](https://doi.org/10.1021/la703434s)

[2] R.V. Listyarini, J. Gamper, T.S. Hofer, **"Storage and Diffusion of Carbon Dioxide in the Metal Organic Framework MOF-5—A Semi-empirical Molecular Dynamics Study,"** *J. Phys. Chem. B*, 2023, 127(43), 9378–9389. DOI: [10.1021/acs.jpcb.3c04155](https://doi.org/10.1021/acs.jpcb.3c04155)

[3] Z. Zhao, Z. Li, Y.S. Lin, **"Adsorption and Diffusion of Carbon Dioxide on Metal–Organic Framework (MOF-5),"** *Ind. Eng. Chem. Res.*, 2009, 48(22), 10015–10020. DOI: [10.1021/ie900665f](https://doi.org/10.1021/ie900665f)

[4] W. Liu, L. Li, J. Liu, L. Ma, C. Fu, Z. Liu, D. Jing, **"Interactions between Absorbed Components of Various CO₂/N₂ Mixtures over MOF-5: Molecular Dynamics Simulation and Density Functional Theory Calculation,"** *Sep. Purif. Technol.*, 2025, 361, 131646. DOI: [10.1016/j.seppur.2025.131646](https://doi.org/10.1016/j.seppur.2025.131646)

[5] S. Cheng *et al.*, **"Multi-scale Design of MOF-based Membrane Separation for CO₂/CH₄ Mixture via Integration of Molecular Simulation, Machine Learning and Process Modeling and Simulation,"** *J. Membr. Sci.*, 2024, 690, 122365. DOI: [10.1016/j.memsci.2023.122365](https://doi.org/10.1016/j.memsci.2023.122365)

[6] D. Saha, Z. Bao, F. Jia, S. Deng, **"Adsorption of CO₂, CH₄, N₂O, and N₂ on MOF-5, MOF-177, and Zeolite 5A,"** *Environ. Sci. Technol.*, 2010, 44(5), 1820–1826. DOI: [10.1021/es9032309](https://doi.org/10.1021/es9032309)

[7] J.-S. Choi, W.-J. Son, J. Kim, W.-S. Ahn, **"Metal–Organic Framework MOF-5 Prepared by Microwave Heating: Factors to be Considered,"** *Microporous Mesoporous Mater.*, 2008, 116(1–3), 727–731. DOI: [10.1016/j.micromeso.2008.04.033](https://doi.org/10.1016/j.micromeso.2008.04.033)

---

## Detailed Analysis

### 1. Is my adsorption capacity reasonable?

**Yes.** Your adsorption capacity of **30.59 mg/g (0.695 mmol/g)** at 298 K and 1 bar is well within the expected range for MOF-5.

- Saha *et al.* (2010) measured CO₂ adsorption on MOF-5 volumetrically up to ~1.05 bar (800 Torr) at 298 K; their isotherm shape suggests low-pressure uptake in the range of **20–40 mg/g** .
- Choi *et al.* (2008) reported a heat of adsorption of **15.8–16.5 kJ/mol** at low pressure (<1 atm), consistent with physisorption in the 0.5–1.5 mmol/g range at ambient pressure .
- The Sarmiento-Pérez GCMC study (cited in the *Theoretical Study of CH₄ and CO₂ Separation by IRMOFs*) reports an adsorption enthalpy of **–13.02 kJ/mol** for CO₂ in IRMOF-1 at 353 K / 1.33 bar .

Your loading of **0.695 mol/kg** corresponds to approximately **1 molecule per primitive cell** or moderate occupancy of the large MOF-5 cavities, which is physically sensible at 1 bar and 298 K.

---

### 2. Is my diffusion coefficient physically reasonable?

**Yes, absolutely.** Your self-diffusion coefficient of **1.616 × 10⁻⁸ m²/s** is physically reasonable and falls within the narrow literature range for MD self-diffusion of CO₂ in MOF-5.

Published MD self-diffusion coefficients for CO₂ in MOF-5/IRMOF-1 at ~298 K:
- **Babarao & Jiang (2008):** 1.4–3.0 × 10⁻⁸ m²/s 
- **Listyarini *et al.* (2023):** 2.09 × 10⁻⁸ m²/s (single CO₂ molecule, 298 K) 

Your value is **within 15% of the lower Babarao & Jiang value** and **within 23% of Listyarini's single-molecule value**. This is excellent agreement given the differences in force fields (UFF vs. DFTB), framework treatment (rigid vs. flexible), and loading.

**Important distinction:** Your MD self-diffusion coefficient should **NOT** be compared to experimental micropore diffusion coefficients:
- Zhao *et al.* (2009) report **8.1–11.5 × 10⁻⁹ cm²/s** (= 8.1–11.5 × 10⁻¹³ m²/s) from gravimetric uptake kinetics .
- Saha *et al.* (2010) report average diffusivities of **~10⁻⁹ m²/s** from uptake kinetics .

These experimental values are **4–5 orders of magnitude smaller** than your MD self-diffusion coefficient because they measure **macroscopic uptake kinetics** dominated by intercrystalline transport, surface barriers, and grain boundaries—not the intrinsic single-molecule self-diffusion measured in MD .

---

### 3. Does my result agree with published literature?

**Yes, with one caveat.**

| Property | Your Result | Literature | Assessment |
|:---|:---|:---|:---|
| **Adsorption capacity** | 30.6 mg/g | 20–40 mg/g (1 bar, 298 K) | ✅ Excellent agreement |
| **MD Self-diffusion** | 1.62 × 10⁻⁸ m²/s | 1.4–3.0 × 10⁻⁸ m²/s | ✅ Excellent agreement |
| **Host–CO₂ interaction** | ~–10.8 kJ/mol | –13.0 kJ/mol (sim.) / –15.1 kJ/mol (exp.) | ⚠️ Somewhat weaker binding |
| **Anisotropy** | 49.4% (Dz << Dx, Dy) | Isotropic expected (cubic Fm-3̄m) | ❌ Red flag |

Your adsorption capacity and self-diffusion coefficient are in excellent agreement with published computational studies. However, your **per-molecule host–CO₂ interaction energy (–10.8 kJ/mol) is ~20% weaker** than the DFTB value (–13.0 kJ/mol) and ~30% weaker than experimental heats of adsorption (–15.1 kJ/mol). This is a known limitation of the **UFF force field with generic charges** for polar adsorbates in MOFs.

---

### 4. Which published paper is MOST similar to my methodology?

**Babarao & Jiang (2008), *Langmuir*** is the most methodologically similar published study [1]. Both use:
- Classical force-field MD (rather than DFT/DFTB)
- **Rigid framework approximation**
- CO₂ self-diffusion in IRMOF-1 via the **Einstein relation**
- Loading-dependent diffusion analysis

**Liu *et al.* (2025), *Sep. Purif. Technol.*** [4] is also highly similar because they employ the **exact same 3×3×2 supercell geometry** and a **GCMC + EMD multiscale workflow** . However, their paper focuses on CO₂/N₂ mixtures and uses DFT for interaction analysis, whereas your study is single-component with classical force fields.

---

### 5. Which published paper is MOST similar to my diffusion coefficient?

**Babarao & Jiang (2008)** [1]. Listyarini *et al.* explicitly cite Babarao & Jiang's CO₂ self-diffusivity in IRMOF-1 as **1.4 × 10⁻⁸ m²/s** . Your value of **1.616 × 10⁻⁸ m²/s** differs by only **~15%** from this literature benchmark. This is remarkably close agreement considering:
- Different force fields (UFF vs. likely DREIDING/UFF hybrid in Babarao & Jiang)
- Your rigid framework vs. their framework treatment
- Different loading conditions

---

### 6. Are there any red flags?

**Yes. Three issues require attention:**

#### 🚩 Red Flag 1: Anisotropy in a Cubic Crystal
MOF-5 crystallizes in the **cubic Fm-3̄m** space group. In a bulk cubic crystal, self-diffusion must be **isotropic** (Dₓ = Dᵧ = D₂). Your reported anisotropy of **49.4%** with D₂ significantly lower than Dₓ/Dᵧ is physically inconsistent with the crystal symmetry.

**Likely causes:**
- **Finite-size effect from the non-cubic 3×3×2 supercell** (55.26 × 55.26 × 36.84 Å). The shorter z-dimension restricts long-wavelength fluctuations and may artificially suppress z-direction diffusion.
- **Insufficient sampling** in the z-direction due to the asymmetric box.
- **Periodic image interactions** in the z-direction.

**Recommended fix:** Re-run with a **cubic supercell** (e.g., 2×2×2 or 3×3×3) or explicitly demonstrate that the anisotropy is a finite-size artifact that vanishes with larger boxes.

#### 🚩 Red Flag 2: Understimated Host–Guest Interaction Energy
Your per-molecule host–CO₂ interaction energy of **–10.8 kJ/mol** is weaker than:
- Listyarini's DFTB value: **–13.0 kJ/mol** 
- Experimental heat of adsorption: **–15.1 to –14.9 kJ/mol** 
- Choi *et al.* experimental Qst: **15.8–16.5 kJ/mol** at low pressure 

This suggests the **UFF force field with standard charges may underestimate CO₂ binding** in MOF-5. Consider validating with DDEC or REPEAT charges, or comparing with the Dubbeldam–Snurr–Vlugt (DSV) MOF force field.

#### 🚩 Red Flag 3: Single Pressure Point
You report results at only **1 bar**. A reviewer will expect at least a **partial isotherm** (e.g., 0.1, 0.5, 1, 5, 10 bar) to demonstrate that your GCMC simulation correctly captures the adsorption behavior across the pressure range and to enable calculation of the isosteric heat of adsorption.

---

### 7. What would a reviewer criticize?

A critical reviewer at a top-tier journal (*J. Phys. Chem. C*, *Langmuir*, *Chem. Eng. J.*) would likely raise the following points:

1. **"Why is diffusion anisotropic in a cubic MOF?"** This is the most serious criticism. The 49.4% anisotropy contradicts the Fm-3̄m symmetry of MOF-5. You must either explain this as a finite-size artifact or correct it.

2. **"Why is the heat of adsorption not reported?"** You provide host–CO₂ interaction energy in Kelvin but do not convert to kJ/mol or compare with experimental isosteric heats. Reviewers expect this conversion and comparison.

3. **"Why is the framework rigid?"** While the rigid approximation is common and justified for MOF-5 at low loading , a reviewer may ask for a justification or a sensitivity test showing that flexibility does not significantly affect Ds at your loading.

4. **"What charge model was used for the framework?"** You mention Ewald summation but do not specify how partial charges were assigned to MOF-5 atoms (e.g., QEq, DDEC, REPEAT, literature values). This is essential for reproducibility.

5. **"Only one pressure point?"** A single data point at 1 bar is insufficient to validate the GCMC methodology. At minimum, provide a 3–5 point isotherm.

6. **"No comparison with experimental uptake kinetics?"** While you correctly avoid comparing MD self-diffusion with experimental micropore diffusion, a reviewer may still ask why you don't discuss the 4–5 order-of-magnitude discrepancy with experimental diffusivities. You should explicitly address this distinction in your manuscript.

7. **"Why 5 ns production?"** While 5 ns is generous, a reviewer might ask for a convergence analysis (e.g., Ds vs. production time) to prove the system is fully sampled.

---

### 8. What are the strengths of my methodology?

Your study has several genuine methodological strengths:

1. **Multiscale GCMC–MD workflow:** The sequential GCMC → NVT → NVE approach is the gold standard for adsorption/diffusion studies and matches the workflow used in high-impact studies like Cheng *et al.* (2024) .

2. **Long production run:** Your **5 ns production** exceeds the 0.53 ns per temperature point used by Listyarini *et al.* (2023) and is comparable to the 5 ns production used in the large-scale MOF membrane screening study .

3. **Block averaging with error bars:** Your block-averaged Ds = (1.700 ± 0.178) × 10⁻⁴ cm²/s provides a proper statistical uncertainty estimate (~10%), which many published studies omit.

4. **Automatic diffusive regime identification:** This is a sophisticated analysis feature that demonstrates rigorous MSD analysis.

5. **COM-based MSD:** Using center-of-mass MSD (rather than atomic MSD) is the correct approach for rigid linear molecules like CO₂.

6. **Ewald summation:** Proper treatment of long-range electrostatics is essential for CO₂ in MOFs and is correctly implemented.

7. **Per-molecule analysis:** Identifying one low-mobility molecule without widespread trapping shows good physical insight.

---

### 9. Can this be considered publication-quality?

**Yes, conditionally.** The core results (adsorption capacity and self-diffusion coefficient) are physically sound and consistent with Q1 literature. The methodology is standard and well-executed. **However, the anisotropy issue must be resolved before submission.**

With the following revisions, this work would be competitive for **Langmuir, J. Phys. Chem. C, or Ind. Eng. Chem. Res.**:
- **Fix or explain the anisotropy** (preferably by re-running with a cubic supercell)
- **Report the isosteric heat of adsorption** (or at least convert interaction energies to kJ/mol)
- **Add 2–4 more pressure points** to generate a partial isotherm
- **Explicitly state the charge assignment method** for MOF-5
- **Add a paragraph distinguishing** MD self-diffusion from experimental uptake kinetics

Without these revisions, a reviewer would likely recommend **major revision**.

---

### 10. Estimate what percentile my methodology falls into compared with published computational MOF diffusion studies.

I estimate your methodology falls into the **~65th–75th percentile** of published computational MOF diffusion studies. Here is the rationale:

| Tier | Percentile | Characteristics | Your Position |
|:---|:---|:---|:---|
| **Top tier** | 85–100th | DFT/DFTB-based MD, flexible frameworks, multiple validation against experiment, systematic force-field comparison, >10 ns production, high-throughput screening | Below this |
| **Upper-middle** | 65–85th | Classical GCMC+MD, proper Ewald, block averaging, error analysis, multiple pressures/loadings, comparison with ≥3 literature sources | **You are here** |
| **Middle** | 35–65th | Classical MD only, no GCMC coupling, short production (<2 ns), no error bars, single pressure, limited literature comparison | Above this |
| **Lower** | <35th | Undersampled MD, incorrect diffusion analysis (e.g., fitting ballistic regime), no electrostatics, comparison only with experimental kinetics | Well above |

**What would push you into the top tier (85th+ percentile):**
- Framework flexibility (NPT or NVT with flexible linkers)
- Systematic force-field comparison (UFF vs. DREIDING vs. BTW-FF)
- Charge sensitivity analysis (QEq vs. DDEC vs. REPEAT)
- Longer production (10–20 ns) with full convergence analysis
- Multiple pressures and temperatures with Arrhenius analysis for activation energy
- Direct comparison with quasi-elastic neutron scattering (QENS) or PFG-NMR data

---

## Reviewer-Style Final Assessment

> **Overall Assessment: GOOD → VERY GOOD (pending revisions)**

| Criterion | Rating | Justification |
|:---|:---|:---|
| **Physical reasonableness of adsorption** | ✅ Very Good | 30.6 mg/g at 1 bar/298 K is consistent with experimental and computational literature. |
| **Physical reasonableness of diffusion** | ✅ Very Good | 1.616 × 10⁻⁸ m²/s is within 15% of the most cited MD benchmark (Babarao & Jiang). |
| **Methodological rigor** | ⚠️ Good | GCMC+MD workflow is standard; 5 ns production is generous; block averaging is proper. Anisotropy and single-pressure data are weaknesses. |
| **Literature comparison** | ⚠️ Good | Comparison with Babarao & Jiang and Listyarini is strong. Missing comparison with experimental Qst and isotherm data. |
| **Reproducibility** | ⚠️ Acceptable | Force field and CO₂ model are stated, but framework charge assignment method is not specified. |
| **Red flags** | ❌ Present | Anisotropy in cubic crystal is a serious issue. Underestimated binding energy is a moderate concern. |

### Recommendation

**Accept with Major Revisions** — The simulation results are fundamentally sound and publication-worthy, but the anisotropy issue must be addressed, the isosteric heat must be reported and compared with literature, and additional pressure points should be added to validate the GCMC methodology. With these revisions, this work would make a solid contribution to the MOF diffusion literature at the level of *J. Phys. Chem. C* or *Langmuir*.I'll conduct a comprehensive literature search to find Q1 papers matching your criteria, then perform the detailed comparison and analysis you've requested.


---

# Comprehensive Literature Comparison: CO₂ Adsorption & Self-Diffusion in MOF-5

## Executive Summary

Your simulation results are **physically reasonable and largely consistent with published Q1 literature** for CO₂ in MOF-5 (IRMOF-1). Your self-diffusion coefficient of **1.616 × 10⁻⁸ m²/s** (1.616 × 10⁻⁴ cm²/s) sits squarely within the narrow band of published MD self-diffusion coefficients (1.4–3.0 × 10⁻⁸ m²/s). Your adsorption capacity of **30.6 mg/g** at 298 K / 1 bar is also consistent with low-pressure GCMC and experimental data. However, there are **methodological red flags**—most notably the pronounced diffusion anisotropy in a cubic crystal—that require attention before journal submission.

---

## Detailed Comparison Table

| Parameter | **Your Work** | Babarao & Jiang (2008) [1] | Listyarini *et al.* (2023) [2] | Zhao *et al.* (2009) [3] | Liu *et al.* (2025) [4] | Cheng *et al.* (2024) [5] | Saha *et al.* (2010) [6] | Choi *et al.* (2008) [7] |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| **Journal** | — | *Langmuir* | *J. Phys. Chem. B* | *Ind. Eng. Chem. Res.* | *Sep. Purif. Technol.* | *J. Membr. Sci.* | *Environ. Sci. Technol.* | *Microporous Mesoporous Mater.* |
| **Quartile** | — | Q1 | Q1 | Q1 | Q1 | Q1 | Q1 | Q1 |
| **Year** | 2026 | 2008 | 2023 | 2009 | 2025 | 2024 | 2010 | 2008 |
| **Force Field** | UFF + TraPPE | Classical FF (DREIDING/UFF-like) | SCC-DFTB/3ob/D3 | — (Experimental) | DFT + MD | UFF / TraPPE | — (Experimental) | — (Experimental) |
| **Framework** | Rigid | Rigid | Flexible (NPT) | Crystalline | Rigid/flexible | Rigid | Crystalline | Crystalline |
| **Temperature** | 298 K | 298–398 K | 248–398 K | 295–331 K | 298 K | 298 K | 298 K | 298 K |
| **Pressure** | 1 bar | Various | 1.013 bar | ≤1 atm | 0.1–0.3 MPa | Various | ≤1.05 bar | 40 bar |
| **CO₂ Model** | TraPPE rigid | — | DFTB-derived | — | — | TraPPE | — | — |
| **Supercell** | 3×3×2 (55.3×55.3×36.8 Å) | — | 1×1×1 (~26.1 Å) | 40–60 μm crystals | 3×3×2 (8.0×8.0×5.6 nm) | ≥24 Å cutoff | — | — |
| **Simulation Method** | GCMC → NVT → NVE MD | MD (NVT/NVE) | NPT MD | Gravimetric | GCMC + EMD + DFT | GCMC + MD | Volumetric/Gravimetric | Gravimetric |
| **Production Length** | 5 ns | — | 0.53 ns / temp | — | — | 5 ns (last) | — | — |
| **Loading (molecules)** | 19 in supercell | Various | 1, 4, 8, 12, 16 | — | Various | Various | — | — |
| **Loading (mol/kg)** | 0.695 | — | — | — | — | — | — | — |
| **Self-Diffusion Coefficient (MD)** | **1.616 × 10⁻⁴ cm²/s** (1.616 × 10⁻⁸ m²/s) | **1.4–3.0 × 10⁻⁸ m²/s** [cited in 2] | **2.09 × 10⁻⁸ m²/s** (1 CO₂) | — | — | Self-diffusivity reported | — | — |
| **Experimental Diffusion** | — | — | — | **8.1–11.5 × 10⁻⁹ cm²/s** (micropore) | — | — | **~10⁻⁹ m²/s** (uptake kinetics) | — |
| **Adsorption Capacity** | 30.59 mg/g (0.695 mmol/g) | — | — | — | — | Consistent with exp. | ~20–35 mg/g (est. at 1 bar) | ~812 mg/g (at 40 bar) |
| **Heat of Adsorption / Interaction Energy** | Host–CO₂: –10.8 kJ/mol (per molecule) | — | –13.0 kJ/mol (int. energy); exp. –15.1 kJ/mol | ~34 kJ/mol (isosteric) | — | Isosteric heat reported | — | 15.8–16.5 kJ/mol (low P) |
| **Activation Energy** | — | 4.05 kJ/mol | 5.47 kJ/mol | 7.61 kJ/mol | — | Reported | — | — |
| **MSD Analysis** | COM-MSD, Einstein, auto diffusive window | MSD, Einstein | MSD, Einstein | — | MSD | MSD | — | — |
| **Anisotropy Reported** | **Dx=1.825, Dy=1.911, Dz=1.112 (49.4%)** | Isotropic (cubic) | Isotropic (cubic) | — | — | — | — | — |
| **Agreement with Your Work** | — | **Excellent** (Ds within 15%) | **Very Good** (Ds within 23%) | N/A (different meas.) | **Methodology match** (same supercell) | Methodology match | Capacity consistent | Qst consistent |

---

### References

[1] R. Babarao, Z. Hu, J. Jiang, **"Diffusion and Separation of CO₂ and CH₄ in Silicalite, C168 Schwarzite, and IRMOF-1: A Comparative Study from Molecular Dynamics Simulation,"** *Langmuir*, 2008, 24(10), 5474–5484. DOI: [10.1021/la703434s](https://doi.org/10.1021/la703434s)

[2] R.V. Listyarini, J. Gamper, T.S. Hofer, **"Storage and Diffusion of Carbon Dioxide in the Metal Organic Framework MOF-5—A Semi-empirical Molecular Dynamics Study,"** *J. Phys. Chem. B*, 2023, 127(43), 9378–9389. DOI: [10.1021/acs.jpcb.3c04155](https://doi.org/10.1021/acs.jpcb.3c04155)

[3] Z. Zhao, Z. Li, Y.S. Lin, **"Adsorption and Diffusion of Carbon Dioxide on Metal–Organic Framework (MOF-5),"** *Ind. Eng. Chem. Res.*, 2009, 48(22), 10015–10020. DOI: [10.1021/ie900665f](https://doi.org/10.1021/ie900665f)

[4] W. Liu, L. Li, J. Liu, L. Ma, C. Fu, Z. Liu, D. Jing, **"Interactions between Absorbed Components of Various CO₂/N₂ Mixtures over MOF-5: Molecular Dynamics Simulation and Density Functional Theory Calculation,"** *Sep. Purif. Technol.*, 2025, 361, 131646. DOI: [10.1016/j.seppur.2025.131646](https://doi.org/10.1016/j.seppur.2025.131646)

[5] S. Cheng *et al.*, **"Multi-scale Design of MOF-based Membrane Separation for CO₂/CH₄ Mixture via Integration of Molecular Simulation, Machine Learning and Process Modeling and Simulation,"** *J. Membr. Sci.*, 2024, 690, 122365. DOI: [10.1016/j.memsci.2023.122365](https://doi.org/10.1016/j.memsci.2023.122365)

[6] D. Saha, Z. Bao, F. Jia, S. Deng, **"Adsorption of CO₂, CH₄, N₂O, and N₂ on MOF-5, MOF-177, and Zeolite 5A,"** *Environ. Sci. Technol.*, 2010, 44(5), 1820–1826. DOI: [10.1021/es9032309](https://doi.org/10.1021/es9032309)

[7] J.-S. Choi, W.-J. Son, J. Kim, W.-S. Ahn, **"Metal–Organic Framework MOF-5 Prepared by Microwave Heating: Factors to be Considered,"** *Microporous Mesoporous Mater.*, 2008, 116(1–3), 727–731. DOI: [10.1016/j.micromeso.2008.04.033](https://doi.org/10.1016/j.micromeso.2008.04.033)

---

## Detailed Analysis

### 1. Is my adsorption capacity reasonable?

**Yes.** Your adsorption capacity of **30.59 mg/g (0.695 mmol/g)** at 298 K and 1 bar is well within the expected range for MOF-5.

- Saha *et al.* (2010) measured CO₂ adsorption on MOF-5 volumetrically up to ~1.05 bar (800 Torr) at 298 K; their isotherm shape suggests low-pressure uptake in the range of **20–40 mg/g** .
- Choi *et al.* (2008) reported a heat of adsorption of **15.8–16.5 kJ/mol** at low pressure (<1 atm), consistent with physisorption in the 0.5–1.5 mmol/g range at ambient pressure .
- The Sarmiento-Pérez GCMC study (cited in the *Theoretical Study of CH₄ and CO₂ Separation by IRMOFs*) reports an adsorption enthalpy of **–13.02 kJ/mol** for CO₂ in IRMOF-1 at 353 K / 1.33 bar .

Your loading of **0.695 mol/kg** corresponds to approximately **1 molecule per primitive cell** or moderate occupancy of the large MOF-5 cavities, which is physically sensible at 1 bar and 298 K.

---

### 2. Is my diffusion coefficient physically reasonable?

**Yes, absolutely.** Your self-diffusion coefficient of **1.616 × 10⁻⁸ m²/s** is physically reasonable and falls within the narrow literature range for MD self-diffusion of CO₂ in MOF-5.

Published MD self-diffusion coefficients for CO₂ in MOF-5/IRMOF-1 at ~298 K:
- **Babarao & Jiang (2008):** 1.4–3.0 × 10⁻⁸ m²/s 
- **Listyarini *et al.* (2023):** 2.09 × 10⁻⁸ m²/s (single CO₂ molecule, 298 K) 

Your value is **within 15% of the lower Babarao & Jiang value** and **within 23% of Listyarini's single-molecule value**. This is excellent agreement given the differences in force fields (UFF vs. DFTB), framework treatment (rigid vs. flexible), and loading.

**Important distinction:** Your MD self-diffusion coefficient should **NOT** be compared to experimental micropore diffusion coefficients:
- Zhao *et al.* (2009) report **8.1–11.5 × 10⁻⁹ cm²/s** (= 8.1–11.5 × 10⁻¹³ m²/s) from gravimetric uptake kinetics .
- Saha *et al.* (2010) report average diffusivities of **~10⁻⁹ m²/s** from uptake kinetics .

These experimental values are **4–5 orders of magnitude smaller** than your MD self-diffusion coefficient because they measure **macroscopic uptake kinetics** dominated by intercrystalline transport, surface barriers, and grain boundaries—not the intrinsic single-molecule self-diffusion measured in MD .

---

### 3. Does my result agree with published literature?

**Yes, with one caveat.**

| Property | Your Result | Literature | Assessment |
|:---|:---|:---|:---|
| **Adsorption capacity** | 30.6 mg/g | 20–40 mg/g (1 bar, 298 K) | ✅ Excellent agreement |
| **MD Self-diffusion** | 1.62 × 10⁻⁸ m²/s | 1.4–3.0 × 10⁻⁸ m²/s | ✅ Excellent agreement |
| **Host–CO₂ interaction** | ~–10.8 kJ/mol | –13.0 kJ/mol (sim.) / –15.1 kJ/mol (exp.) | ⚠️ Somewhat weaker binding |
| **Anisotropy** | 49.4% (Dz << Dx, Dy) | Isotropic expected (cubic Fm-3̄m) | ❌ Red flag |

Your adsorption capacity and self-diffusion coefficient are in excellent agreement with published computational studies. However, your **per-molecule host–CO₂ interaction energy (–10.8 kJ/mol) is ~20% weaker** than the DFTB value (–13.0 kJ/mol) and ~30% weaker than experimental heats of adsorption (–15.1 kJ/mol). This is a known limitation of the **UFF force field with generic charges** for polar adsorbates in MOFs.

---

### 4. Which published paper is MOST similar to my methodology?

**Babarao & Jiang (2008), *Langmuir*** is the most methodologically similar published study [1]. Both use:
- Classical force-field MD (rather than DFT/DFTB)
- **Rigid framework approximation**
- CO₂ self-diffusion in IRMOF-1 via the **Einstein relation**
- Loading-dependent diffusion analysis

**Liu *et al.* (2025), *Sep. Purif. Technol.*** [4] is also highly similar because they employ the **exact same 3×3×2 supercell geometry** and a **GCMC + EMD multiscale workflow** . However, their paper focuses on CO₂/N₂ mixtures and uses DFT for interaction analysis, whereas your study is single-component with classical force fields.

---

### 5. Which published paper is MOST similar to my diffusion coefficient?

**Babarao & Jiang (2008)** [1]. Listyarini *et al.* explicitly cite Babarao & Jiang's CO₂ self-diffusivity in IRMOF-1 as **1.4 × 10⁻⁸ m²/s** . Your value of **1.616 × 10⁻⁸ m²/s** differs by only **~15%** from this literature benchmark. This is remarkably close agreement considering:
- Different force fields (UFF vs. likely DREIDING/UFF hybrid in Babarao & Jiang)
- Your rigid framework vs. their framework treatment
- Different loading conditions

---

### 6. Are there any red flags?

**Yes. Three issues require attention:**

#### 🚩 Red Flag 1: Anisotropy in a Cubic Crystal
MOF-5 crystallizes in the **cubic Fm-3̄m** space group. In a bulk cubic crystal, self-diffusion must be **isotropic** (Dₓ = Dᵧ = D₂). Your reported anisotropy of **49.4%** with D₂ significantly lower than Dₓ/Dᵧ is physically inconsistent with the crystal symmetry.

**Likely causes:**
- **Finite-size effect from the non-cubic 3×3×2 supercell** (55.26 × 55.26 × 36.84 Å). The shorter z-dimension restricts long-wavelength fluctuations and may artificially suppress z-direction diffusion.
- **Insufficient sampling** in the z-direction due to the asymmetric box.
- **Periodic image interactions** in the z-direction.

**Recommended fix:** Re-run with a **cubic supercell** (e.g., 2×2×2 or 3×3×3) or explicitly demonstrate that the anisotropy is a finite-size artifact that vanishes with larger boxes.

#### 🚩 Red Flag 2: Understimated Host–Guest Interaction Energy
Your per-molecule host–CO₂ interaction energy of **–10.8 kJ/mol** is weaker than:
- Listyarini's DFTB value: **–13.0 kJ/mol** 
- Experimental heat of adsorption: **–15.1 to –14.9 kJ/mol** 
- Choi *et al.* experimental Qst: **15.8–16.5 kJ/mol** at low pressure 

This suggests the **UFF force field with standard charges may underestimate CO₂ binding** in MOF-5. Consider validating with DDEC or REPEAT charges, or comparing with the Dubbeldam–Snurr–Vlugt (DSV) MOF force field.

#### 🚩 Red Flag 3: Single Pressure Point
You report results at only **1 bar**. A reviewer will expect at least a **partial isotherm** (e.g., 0.1, 0.5, 1, 5, 10 bar) to demonstrate that your GCMC simulation correctly captures the adsorption behavior across the pressure range and to enable calculation of the isosteric heat of adsorption.

---

### 7. What would a reviewer criticize?

A critical reviewer at a top-tier journal (*J. Phys. Chem. C*, *Langmuir*, *Chem. Eng. J.*) would likely raise the following points:

1. **"Why is diffusion anisotropic in a cubic MOF?"** This is the most serious criticism. The 49.4% anisotropy contradicts the Fm-3̄m symmetry of MOF-5. You must either explain this as a finite-size artifact or correct it.

2. **"Why is the heat of adsorption not reported?"** You provide host–CO₂ interaction energy in Kelvin but do not convert to kJ/mol or compare with experimental isosteric heats. Reviewers expect this conversion and comparison.

3. **"Why is the framework rigid?"** While the rigid approximation is common and justified for MOF-5 at low loading , a reviewer may ask for a justification or a sensitivity test showing that flexibility does not significantly affect Ds at your loading.

4. **"What charge model was used for the framework?"** You mention Ewald summation but do not specify how partial charges were assigned to MOF-5 atoms (e.g., QEq, DDEC, REPEAT, literature values). This is essential for reproducibility.

5. **"Only one pressure point?"** A single data point at 1 bar is insufficient to validate the GCMC methodology. At minimum, provide a 3–5 point isotherm.

6. **"No comparison with experimental uptake kinetics?"** While you correctly avoid comparing MD self-diffusion with experimental micropore diffusion, a reviewer may still ask why you don't discuss the 4–5 order-of-magnitude discrepancy with experimental diffusivities. You should explicitly address this distinction in your manuscript.

7. **"Why 5 ns production?"** While 5 ns is generous, a reviewer might ask for a convergence analysis (e.g., Ds vs. production time) to prove the system is fully sampled.

---

### 8. What are the strengths of my methodology?

Your study has several genuine methodological strengths:

1. **Multiscale GCMC–MD workflow:** The sequential GCMC → NVT → NVE approach is the gold standard for adsorption/diffusion studies and matches the workflow used in high-impact studies like Cheng *et al.* (2024) .

2. **Long production run:** Your **5 ns production** exceeds the 0.53 ns per temperature point used by Listyarini *et al.* (2023) and is comparable to the 5 ns production used in the large-scale MOF membrane screening study .

3. **Block averaging with error bars:** Your block-averaged Ds = (1.700 ± 0.178) × 10⁻⁴ cm²/s provides a proper statistical uncertainty estimate (~10%), which many published studies omit.

4. **Automatic diffusive regime identification:** This is a sophisticated analysis feature that demonstrates rigorous MSD analysis.

5. **COM-based MSD:** Using center-of-mass MSD (rather than atomic MSD) is the correct approach for rigid linear molecules like CO₂.

6. **Ewald summation:** Proper treatment of long-range electrostatics is essential for CO₂ in MOFs and is correctly implemented.

7. **Per-molecule analysis:** Identifying one low-mobility molecule without widespread trapping shows good physical insight.

---

### 9. Can this be considered publication-quality?

**Yes, conditionally.** The core results (adsorption capacity and self-diffusion coefficient) are physically sound and consistent with Q1 literature. The methodology is standard and well-executed. **However, the anisotropy issue must be resolved before submission.**

With the following revisions, this work would be competitive for **Langmuir, J. Phys. Chem. C, or Ind. Eng. Chem. Res.**:
- **Fix or explain the anisotropy** (preferably by re-running with a cubic supercell)
- **Report the isosteric heat of adsorption** (or at least convert interaction energies to kJ/mol)
- **Add 2–4 more pressure points** to generate a partial isotherm
- **Explicitly state the charge assignment method** for MOF-5
- **Add a paragraph distinguishing** MD self-diffusion from experimental uptake kinetics

Without these revisions, a reviewer would likely recommend **major revision**.

---

### 10. Estimate what percentile my methodology falls into compared with published computational MOF diffusion studies.

I estimate your methodology falls into the **~65th–75th percentile** of published computational MOF diffusion studies. Here is the rationale:

| Tier | Percentile | Characteristics | Your Position |
|:---|:---|:---|:---|
| **Top tier** | 85–100th | DFT/DFTB-based MD, flexible frameworks, multiple validation against experiment, systematic force-field comparison, >10 ns production, high-throughput screening | Below this |
| **Upper-middle** | 65–85th | Classical GCMC+MD, proper Ewald, block averaging, error analysis, multiple pressures/loadings, comparison with ≥3 literature sources | **You are here** |
| **Middle** | 35–65th | Classical MD only, no GCMC coupling, short production (<2 ns), no error bars, single pressure, limited literature comparison | Above this |
| **Lower** | <35th | Undersampled MD, incorrect diffusion analysis (e.g., fitting ballistic regime), no electrostatics, comparison only with experimental kinetics | Well above |

**What would push you into the top tier (85th+ percentile):**
- Framework flexibility (NPT or NVT with flexible linkers)
- Systematic force-field comparison (UFF vs. DREIDING vs. BTW-FF)
- Charge sensitivity analysis (QEq vs. DDEC vs. REPEAT)
- Longer production (10–20 ns) with full convergence analysis
- Multiple pressures and temperatures with Arrhenius analysis for activation energy
- Direct comparison with quasi-elastic neutron scattering (QENS) or PFG-NMR data

---

## Reviewer-Style Final Assessment

> **Overall Assessment: GOOD → VERY GOOD (pending revisions)**

| Criterion | Rating | Justification |
|:---|:---|:---|
| **Physical reasonableness of adsorption** | ✅ Very Good | 30.6 mg/g at 1 bar/298 K is consistent with experimental and computational literature. |
| **Physical reasonableness of diffusion** | ✅ Very Good | 1.616 × 10⁻⁸ m²/s is within 15% of the most cited MD benchmark (Babarao & Jiang). |
| **Methodological rigor** | ⚠️ Good | GCMC+MD workflow is standard; 5 ns production is generous; block averaging is proper. Anisotropy and single-pressure data are weaknesses. |
| **Literature comparison** | ⚠️ Good | Comparison with Babarao & Jiang and Listyarini is strong. Missing comparison with experimental Qst and isotherm data. |
| **Reproducibility** | ⚠️ Acceptable | Force field and CO₂ model are stated, but framework charge assignment method is not specified. |
| **Red flags** | ❌ Present | Anisotropy in cubic crystal is a serious issue. Underestimated binding energy is a moderate concern. |

### Recommendation

**Accept with Major Revisions** — The simulation results are fundamentally sound and publication-worthy, but the anisotropy issue must be addressed, the isosteric heat must be reported and compared with literature, and additional pressure points should be added to validate the GCMC methodology. With these revisions, this work would make a solid contribution to the MOF diffusion literature at the level of *J. Phys. Chem. C* or *Langmuir*.
---

## Detailed Comparison Table

| Parameter | **Your Work** | Babarao & Jiang (2008) [1] | Listyarini *et al.* (2023) [2] | Zhao *et al.* (2009) [3] | Liu *et al.* (2025) [4] | Cheng *et al.* (2024) [5] | Saha *et al.* (2010) [6] | Choi *et al.* (2008) [7] |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| **Journal** | — | *Langmuir* | *J. Phys. Chem. B* | *Ind. Eng. Chem. Res.* | *Sep. Purif. Technol.* | *J. Membr. Sci.* | *Environ. Sci. Technol.* | *Microporous Mesoporous Mater.* |
| **Quartile** | — | Q1 | Q1 | Q1 | Q1 | Q1 | Q1 | Q1 |
| **Year** | 2026 | 2008 | 2023 | 2009 | 2025 | 2024 | 2010 | 2008 |
| **Force Field** | UFF + TraPPE | Classical FF (DREIDING/UFF-like) | SCC-DFTB/3ob/D3 | — (Experimental) | DFT + MD | UFF / TraPPE | — (Experimental) | — (Experimental) |
| **Framework** | Rigid | Rigid | Flexible (NPT) | Crystalline | Rigid/flexible | Rigid | Crystalline | Crystalline |
| **Temperature** | 298 K | 298–398 K | 248–398 K | 295–331 K | 298 K | 298 K | 298 K | 298 K |
| **Pressure** | 1 bar | Various | 1.013 bar | ≤1 atm | 0.1–0.3 MPa | Various | ≤1.05 bar | 40 bar |
| **CO₂ Model** | TraPPE rigid | — | DFTB-derived | — | — | TraPPE | — | — |
| **Supercell** | 3×3×2 (55.3×55.3×36.8 Å) | — | 1×1×1 (~26.1 Å) | 40–60 μm crystals | 3×3×2 (8.0×8.0×5.6 nm) | ≥24 Å cutoff | — | — |
| **Simulation Method** | GCMC → NVT → NVE MD | MD (NVT/NVE) | NPT MD | Gravimetric | GCMC + EMD + DFT | GCMC + MD | Volumetric/Gravimetric | Gravimetric |
| **Production Length** | 5 ns | — | 0.53 ns / temp | — | — | 5 ns (last) | — | — |
| **Loading (molecules)** | 19 in supercell | Various | 1, 4, 8, 12, 16 | — | Various | Various | — | — |
| **Loading (mol/kg)** | 0.695 | — | — | — | — | — | — | — |
| **Self-Diffusion Coefficient (MD)** | **1.616 × 10⁻⁴ cm²/s** (1.616 × 10⁻⁸ m²/s) | **1.4–3.0 × 10⁻⁸ m²/s** [cited in 2] | **2.09 × 10⁻⁸ m²/s** (1 CO₂) | — | — | Self-diffusivity reported | — | — |
| **Experimental Diffusion** | — | — | — | **8.1–11.5 × 10⁻⁹ cm²/s** (micropore) | — | — | **~10⁻⁹ m²/s** (uptake kinetics) | — |
| **Adsorption Capacity** | 30.59 mg/g (0.695 mmol/g) | — | — | — | — | Consistent with exp. | ~20–35 mg/g (est. at 1 bar) | ~812 mg/g (at 40 bar) |
| **Heat of Adsorption / Interaction Energy** | Host–CO₂: –10.8 kJ/mol (per molecule) | — | –13.0 kJ/mol (int. energy); exp. –15.1 kJ/mol | ~34 kJ/mol (isosteric) | — | Isosteric heat reported | — | 15.8–16.5 kJ/mol (low P) |
| **Activation Energy** | — | 4.05 kJ/mol | 5.47 kJ/mol | 7.61 kJ/mol | — | Reported | — | — |
| **MSD Analysis** | COM-MSD, Einstein, auto diffusive window | MSD, Einstein | MSD, Einstein | — | MSD | MSD | — | — |
| **Anisotropy Reported** | **Dx=1.825, Dy=1.911, Dz=1.112 (49.4%)** | Isotropic (cubic) | Isotropic (cubic) | — | — | — | — | — |
| **Agreement with Your Work** | — | **Excellent** (Ds within 15%) | **Very Good** (Ds within 23%) | N/A (different meas.) | **Methodology match** (same supercell) | Methodology match | Capacity consistent | Qst consistent |

---

### References

[1] R. Babarao, Z. Hu, J. Jiang, **"Diffusion and Separation of CO₂ and CH₄ in Silicalite, C168 Schwarzite, and IRMOF-1: A Comparative Study from Molecular Dynamics Simulation,"** *Langmuir*, 2008, 24(10), 5474–5484. DOI: [10.1021/la703434s](https://doi.org/10.1021/la703434s)

[2] R.V. Listyarini, J. Gamper, T.S. Hofer, **"Storage and Diffusion of Carbon Dioxide in the Metal Organic Framework MOF-5—A Semi-empirical Molecular Dynamics Study,"** *J. Phys. Chem. B*, 2023, 127(43), 9378–9389. DOI: [10.1021/acs.jpcb.3c04155](https://doi.org/10.1021/acs.jpcb.3c04155)

[3] Z. Zhao, Z. Li, Y.S. Lin, **"Adsorption and Diffusion of Carbon Dioxide on Metal–Organic Framework (MOF-5),"** *Ind. Eng. Chem. Res.*, 2009, 48(22), 10015–10020. DOI: [10.1021/ie900665f](https://doi.org/10.1021/ie900665f)

[4] W. Liu, L. Li, J. Liu, L. Ma, C. Fu, Z. Liu, D. Jing, **"Interactions between Absorbed Components of Various CO₂/N₂ Mixtures over MOF-5: Molecular Dynamics Simulation and Density Functional Theory Calculation,"** *Sep. Purif. Technol.*, 2025, 361, 131646. DOI: [10.1016/j.seppur.2025.131646](https://doi.org/10.1016/j.seppur.2025.131646)

[5] S. Cheng *et al.*, **"Multi-scale Design of MOF-based Membrane Separation for CO₂/CH₄ Mixture via Integration of Molecular Simulation, Machine Learning and Process Modeling and Simulation,"** *J. Membr. Sci.*, 2024, 690, 122365. DOI: [10.1016/j.memsci.2023.122365](https://doi.org/10.1016/j.memsci.2023.122365)

[6] D. Saha, Z. Bao, F. Jia, S. Deng, **"Adsorption of CO₂, CH₄, N₂O, and N₂ on MOF-5, MOF-177, and Zeolite 5A,"** *Environ. Sci. Technol.*, 2010, 44(5), 1820–1826. DOI: [10.1021/es9032309](https://doi.org/10.1021/es9032309)

[7] J.-S. Choi, W.-J. Son, J. Kim, W.-S. Ahn, **"Metal–Organic Framework MOF-5 Prepared by Microwave Heating: Factors to be Considered,"** *Microporous Mesoporous Mater.*, 2008, 116(1–3), 727–731. DOI: [10.1016/j.micromeso.2008.04.033](https://doi.org/10.1016/j.micromeso.2008.04.033)

---

## Detailed Analysis

### 1. Is my adsorption capacity reasonable?

**Yes.** Your adsorption capacity of **30.59 mg/g (0.695 mmol/g)** at 298 K and 1 bar is well within the expected range for MOF-5.

- Saha *et al.* (2010) measured CO₂ adsorption on MOF-5 volumetrically up to ~1.05 bar (800 Torr) at 298 K; their isotherm shape suggests low-pressure uptake in the range of **20–40 mg/g**.
- Choi *et al.* (2008) reported a heat of adsorption of **15.8–16.5 kJ/mol** at low pressure (<1 atm), consistent with physisorption in the 0.5–1.5 mmol/g range at ambient pressure.
- The Sarmiento-Pérez GCMC study (cited in the *Theoretical Study of CH₄ and CO₂ Separation by IRMOFs*) reports an adsorption enthalpy of **–13.02 kJ/mol** for CO₂ in IRMOF-1 at 353 K / 1.33 bar.

Your loading of **0.695 mol/kg** corresponds to approximately **1 molecule per primitive cell** or moderate occupancy of the large MOF-5 cavities, which is physically sensible at 1 bar and 298 K.

---

### 2. Is my diffusion coefficient physically reasonable?

**Yes, absolutely.** Your self-diffusion coefficient of **1.616 × 10⁻⁸ m²/s** is physically reasonable and falls within the narrow literature range for MD self-diffusion of CO₂ in MOF-5.

Published MD self-diffusion coefficients for CO₂ in MOF-5/IRMOF-1 at ~298 K:
- **Babarao & Jiang (2008):** 1.4–3.0 × 10⁻⁸ m²/s 
- **Listyarini *et al.* (2023):** 2.09 × 10⁻⁸ m²/s (single CO₂ molecule, 298 K) 

Your value is **within 15% of the lower Babarao & Jiang value** and **within 23% of Listyarini's single-molecule value**. This is excellent agreement given the differences in force fields (UFF vs. DFTB), framework treatment (rigid vs. flexible), and loading.

**Important distinction:** Your MD self-diffusion coefficient should **NOT** be compared to experimental micropore diffusion coefficients:
- Zhao *et al.* (2009) report **8.1–11.5 × 10⁻⁹ cm²/s** (= 8.1–11.5 × 10⁻¹³ m²/s) from gravimetric uptake kinetics.
- Saha *et al.* (2010) report average diffusivities of **~10⁻⁹ m²/s** from uptake kinetics.

These experimental values are **4–5 orders of magnitude smaller** than your MD self-diffusion coefficient because they measure **macroscopic uptake kinetics** dominated by intercrystalline transport, surface barriers, and grain boundaries—not the intrinsic single-molecule self-diffusion measured in MD.

---

### 3. Does my result agree with published literature?

**Yes, with one caveat.**

| Property | Your Result | Literature | Assessment |
|:---|:---|:---|:---|
| **Adsorption capacity** | 30.6 mg/g | 20–40 mg/g (1 bar, 298 K) | ✅ Excellent agreement |
| **MD Self-diffusion** | 1.62 × 10⁻⁸ m²/s | 1.4–3.0 × 10⁻⁸ m²/s | ✅ Excellent agreement |
| **Host–CO₂ interaction** | ~–10.8 kJ/mol | –13.0 kJ/mol (sim.) / –15.1 kJ/mol (exp.) | ⚠️ Somewhat weaker binding |
| **Anisotropy** | 49.4% (Dz << Dx, Dy) | Isotropic expected (cubic Fm-3̄m) | ❌ Red flag |

Your adsorption capacity and self-diffusion coefficient are in excellent agreement with published computational studies. However, your **per-molecule host–CO₂ interaction energy (–10.8 kJ/mol) is ~20% weaker** than the DFTB value (–13.0 kJ/mol) and ~30% weaker than experimental heats of adsorption (–15.1 kJ/mol). This is a known limitation of the **UFF force field with generic charges** for polar adsorbates in MOFs.

---

### 4. Which published paper is MOST similar to my methodology?

**Babarao & Jiang (2008), *Langmuir*** is the most methodologically similar published study [1]. Both use:
- Classical force-field MD (rather than DFT/DFTB)
- **Rigid framework approximation**
- CO₂ self-diffusion in IRMOF-1 via the **Einstein relation**
- Loading-dependent diffusion analysis

**Liu *et al.* (2025), *Sep. Purif. Technol.*** [4] is also highly similar because they employ the **exact same 3×3×2 supercell geometry** and a **GCMC + EMD multiscale workflow**. However, their paper focuses on CO₂/N₂ mixtures and uses DFT for interaction analysis, whereas your study is single-component with classical force fields.

---

### 5. Which published paper is MOST similar to my diffusion coefficient?

**Babarao & Jiang (2008)** [1]. Listyarini *et al.* explicitly cite Babarao & Jiang's CO₂ self-diffusivity in IRMOF-1 as **1.4 × 10⁻⁸ m²/s**. Your value of **1.616 × 10⁻⁸ m²/s** differs by only **~15%** from this literature benchmark. This is remarkably close agreement considering:
- Different force fields (UFF vs. likely DREIDING/UFF hybrid in Babarao & Jiang)
- Your rigid framework vs. their framework treatment
- Different loading conditions

---

### 6. Are there any red flags?

**Yes. Three issues require attention:**

#### 🚩 Red Flag 1: Anisotropy in a Cubic Crystal
MOF-5 crystallizes in the **cubic Fm-3̄m** space group. In a bulk cubic crystal, self-diffusion must be **isotropic** (Dₓ = Dᵧ = D₂). Your reported anisotropy of **49.4%** with D₂ significantly lower than Dₓ/Dᵧ is physically inconsistent with the crystal symmetry.

**Likely causes:**
- **Finite-size effect from the non-cubic 3×3×2 supercell** (55.26 × 55.26 × 36.84 Å). The shorter z-dimension restricts long-wavelength fluctuations and may artificially suppress z-direction diffusion.
- **Insufficient sampling** in the z-direction due to the asymmetric box.
- **Periodic image interactions** in the z-direction.

**Recommended fix:** Re-run with a **cubic supercell** (e.g., 2×2×2 or 3×3×3) or explicitly demonstrate that the anisotropy is a finite-size artifact that vanishes with larger boxes.

#### 🚩 Red Flag 2: Understimated Host–Guest Interaction Energy
Your per-molecule host–CO₂ interaction energy of **–10.8 kJ/mol** is weaker than:
- Listyarini's DFTB value: **–13.0 kJ/mol** 
- Experimental heat of adsorption: **–15.1 to –14.9 kJ/mol** 
- Choi *et al.* experimental Qst: **15.8–16.5 kJ/mol** at low pressure 

This suggests the **UFF force field with standard charges may underestimate CO₂ binding** in MOF-5. Consider validating with DDEC or REPEAT charges, or comparing with the Dubbeldam–Snurr–Vlugt (DSV) MOF force field.

#### 🚩 Red Flag 3: Single Pressure Point
You report results at only **1 bar**. A reviewer will expect at least a **partial isotherm** (e.g., 0.1, 0.5, 1, 5, 10 bar) to demonstrate that your GCMC simulation correctly captures the adsorption behavior across the pressure range and to enable calculation of the isosteric heat of adsorption.

---

### 7. What would a reviewer criticize?

A critical reviewer at a top-tier journal (*J. Phys. Chem. C*, *Langmuir*, *Chem. Eng. J.*) would likely raise the following points:

1. **"Why is diffusion anisotropic in a cubic MOF?"** This is the most serious criticism. The 49.4% anisotropy contradicts the Fm-3̄m symmetry of MOF-5. You must either explain this as a finite-size artifact or correct it.

2. **"Why is the heat of adsorption not reported?"** You provide host–CO₂ interaction energy in Kelvin but do not convert to kJ/mol or compare with experimental isosteric heats. Reviewers expect this conversion and comparison.

3. **"Why is the framework rigid?"** While the rigid approximation is common and justified for MOF-5 at low loading, a reviewer may ask for a justification or a sensitivity test showing that flexibility does not significantly affect Ds at your loading.

4. **"What charge model was used for the framework?"** You mention Ewald summation but do not specify how partial charges were assigned to MOF-5 atoms (e.g., QEq, DDEC, REPEAT, literature values). This is essential for reproducibility.

5. **"Only one pressure point?"** A single data point at 1 bar is insufficient to validate the GCMC methodology. At minimum, provide a 3–5 point isotherm.

6. **"No comparison with experimental uptake kinetics?"** While you correctly avoid comparing MD self-diffusion with experimental micropore diffusion, a reviewer may still ask why you don't discuss the 4–5 order-of-magnitude discrepancy with experimental diffusivities. You should explicitly address this distinction in your manuscript.

7. **"Why 5 ns production?"** While 5 ns is generous, a reviewer might ask for a convergence analysis (e.g., Ds vs. production time) to prove the system is fully sampled.

---

### 8. What are the strengths of my methodology?

Your study has several genuine methodological strengths:

1. **Multiscale GCMC–MD workflow:** The sequential GCMC → NVT → NVE approach is the gold standard for adsorption/diffusion studies and matches the workflow used in high-impact studies like Cheng *et al.* (2024).

2. **Long production run:** Your **5 ns production** exceeds the 0.53 ns per temperature point used by Listyarini *et al.* (2023) and is comparable to the 5 ns production used in the large-scale MOF membrane screening study.

3. **Block averaging with error bars:** Your block-averaged Ds = (1.700 ± 0.178) × 10⁻⁴ cm²/s provides a proper statistical uncertainty estimate (~10%), which many published studies omit.

4. **Automatic diffusive regime identification:** This is a sophisticated analysis feature that demonstrates rigorous MSD analysis.

5. **COM-based MSD:** Using center-of-mass MSD (rather than atomic MSD) is the correct approach for rigid linear molecules like CO₂.

6. **Ewald summation:** Proper treatment of long-range electrostatics is essential for CO₂ in MOFs and is correctly implemented.

7. **Per-molecule analysis:** Identifying one low-mobility molecule without widespread trapping shows good physical insight.

---

### 9. Can this be considered publication-quality?

**Yes, conditionally.** The core results (adsorption capacity and self-diffusion coefficient) are physically sound and consistent with Q1 literature. The methodology is standard and well-executed. **However, the anisotropy issue must be resolved before submission.**

With the following revisions, this work would be competitive for **Langmuir, J. Phys. Chem. C, or Ind. Eng. Chem. Res.**:
- **Fix or explain the anisotropy** (preferably by re-running with a cubic supercell)
- **Report the isosteric heat of adsorption** (or at least convert interaction energies to kJ/mol)
- **Add 2–4 more pressure points** to generate a partial isotherm
- **Explicitly state the charge assignment method** for MOF-5
- **Add a paragraph distinguishing** MD self-diffusion from experimental uptake kinetics

Without these revisions, a reviewer would likely recommend **major revision**.

---

### 10. Estimate what percentile my methodology falls into compared with published computational MOF diffusion studies.

I estimate your methodology falls into the **~65th–75th percentile** of published computational MOF diffusion studies. Here is the rationale:

| Tier | Percentile | Characteristics | Your Position |
|:---|:---|:---|:---|
| **Top tier** | 85–100th | DFT/DFTB-based MD, flexible frameworks, multiple validation against experiment, systematic force-field comparison, >10 ns production, high-throughput screening | Below this |
| **Upper-middle** | 65–85th | Classical GCMC+MD, proper Ewald, block averaging, error analysis, multiple pressures/loadings, comparison with ≥3 literature sources | **You are here** |
| **Middle** | 35–65th | Classical MD only, no GCMC coupling, short production (<2 ns), no error bars, single pressure, limited literature comparison | Above this |
| **Lower** | <35th | Undersampled MD, incorrect diffusion analysis (e.g., fitting ballistic regime), no electrostatics, comparison only with experimental kinetics | Well above |

**What would push you into the top tier (85th+ percentile):**
- Framework flexibility (NPT or NVT with flexible linkers)
- Systematic force-field comparison (UFF vs. DREIDING vs. BTW-FF)
- Charge sensitivity analysis (QEq vs. DDEC vs. REPEAT)
- Longer production (10–20 ns) with full convergence analysis
- Multiple pressures and temperatures with Arrhenius analysis for activation energy
- Direct comparison with quasi-elastic neutron scattering (QENS) or PFG-NMR data

---

## Reviewer-Style Final Assessment

> **Overall Assessment: GOOD → VERY GOOD (pending revisions)**

| Criterion | Rating | Justification |
|:---|:---|:---|
| **Physical reasonableness of adsorption** | ✅ Very Good | 30.6 mg/g at 1 bar/298 K is consistent with experimental and computational literature. |
| **Physical reasonableness of diffusion** | ✅ Very Good | 1.616 × 10⁻⁸ m²/s is within 15% of the most cited MD benchmark (Babarao & Jiang). |
| **Methodological rigor** | ⚠️ Good | GCMC+MD workflow is standard; 5 ns production is generous; block averaging is proper. Anisotropy and single-pressure data are weaknesses. |
| **Literature comparison** | ⚠️ Good | Comparison with Babarao & Jiang and Listyarini is strong. Missing comparison with experimental Qst and isotherm data. |
| **Reproducibility** | ⚠️ Acceptable | Force field and CO₂ model are stated, but framework charge assignment method is not specified. |
| **Red flags** | ❌ Present | Anisotropy in cubic crystal is a serious issue. Underestimated binding energy is a moderate concern. |

### Recommendation

**Accept with Major Revisions** — The simulation results are fundamentally sound and publication-worthy, but the anisotropy issue must be addressed, the isosteric heat must be reported and compared with literature, and additional pressure points should be added to validate the GCMC methodology. With these revisions, this work would make a solid contribution to the MOF diffusion literature at the level of *J. Phys. Chem. C* or *Langmuir*.