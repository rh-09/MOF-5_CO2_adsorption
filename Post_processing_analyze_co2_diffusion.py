#!/usr/bin/env python3
"""
CO2 self-diffusion analysis in MOF-5, from a LAMMPS production run.

Reads a LAMMPS custom dump file (columns: id mol type xu yu zu), builds
mass-weighted center-of-mass trajectories for each CO2 molecule, computes
the mean-squared displacement (MSD) via the Einstein relation using an
FFT-accelerated algorithm, automatically identifies the diffusive
(Fickian) regime from the local log-log slope, fits the self-diffusion
coefficient with block-averaged uncertainty, checks for anisotropy and
per-molecule trapping, and writes out plots + a text summary.

Usage:
    python analyze_co2_diffusion.py

Edit the SETTINGS block below to match your actual input script before
running. Requires: numpy, matplotlib.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

# =====================================================================
# SETTINGS -- must match your Stage 3 (production) LAMMPS input script
# =====================================================================
DUMP_FILE      = "prod_unwrapped.dump"
OUT_DIR        = "diffusion_analysis"
TIMESTEP_FS    = 1.0                       # fs per MD step
DUMP_EVERY     = 1000                      # steps between dump frames
N_BLOCKS       = 4                         # blocks for uncertainty estimate
SLOPE_TOL      = 0.15                      # local log-log slope must be within 1 +/- this
SKIP_FRACTION  = 0.10                      # ignore this fraction of run as early/ballistic transient
MASS_MAP       = {6: 12.0107, 7: 15.9994}  # amu: {C_co2: mass, O_co2: mass}
# =====================================================================


def parse_dump(path):
    """Parse a LAMMPS custom dump file with columns: id mol type xu yu zu."""
    timesteps, frames = [], []
    with open(path) as f:
        while True:
            line = f.readline()
            if not line:
                break
            if line.startswith("ITEM: TIMESTEP"):
                timesteps.append(int(f.readline().strip()))
            elif line.startswith("ITEM: NUMBER OF ATOMS"):
                natoms = int(f.readline().strip())
            elif line.startswith("ITEM: BOX BOUNDS"):
                for _ in range(3):
                    f.readline()
            elif line.startswith("ITEM: ATOMS"):
                cols = line.split()[2:]
                idx = {c: i for i, c in enumerate(cols)}
                required = ["id", "mol", "type", "xu", "yu", "zu"]
                missing = [c for c in required if c not in idx]
                if missing:
                    raise ValueError(
                        f"Dump file is missing required columns {missing}. "
                        f"Found columns: {cols}. Did you dump 'xu yu zu' "
                        f"(unwrapped) rather than 'x y z'?"
                    )
                data = np.zeros((natoms, 6))
                for i in range(natoms):
                    parts = f.readline().split()
                    data[i, 0] = float(parts[idx["id"]])
                    data[i, 1] = float(parts[idx["mol"]])
                    data[i, 2] = float(parts[idx["type"]])
                    data[i, 3] = float(parts[idx["xu"]])
                    data[i, 4] = float(parts[idx["yu"]])
                    data[i, 5] = float(parts[idx["zu"]])
                frames.append(data)
    if not frames:
        raise ValueError(f"No frames parsed from {path} -- is the file empty or malformed?")
    return timesteps, frames


def compute_com_trajectories(timesteps, frames, mass_map, timestep_fs):
    """Mass-weighted center-of-mass trajectory per CO2 molecule, per frame.

    Filters to atom types present in mass_map (i.e. CO2 atom types) first --
    the dump may contain the whole system (framework + CO2) rather than just
    the CO2 group, in which case framework atoms/molecule-ids must be
    dropped before building per-molecule center-of-mass trajectories.
    """
    n_frames = len(frames)
    co2_types = set(mass_map.keys())
    frame0_co2 = frames[0][np.isin(frames[0][:, 2].astype(int), list(co2_types))]
    mol_ids = np.unique(frame0_co2[:, 1]).astype(int)
    n_mol = len(mol_ids)
    if n_mol == 0:
        raise ValueError(
            "No atoms matching MASS_MAP types found in the dump -- check "
            "that MASS_MAP keys match your CO2 atom type numbers."
        )
    com = np.zeros((n_frames, n_mol, 3))
    for fi, data in enumerate(frames):
        data = data[np.isin(data[:, 2].astype(int), list(co2_types))]
        order = np.argsort(data[:, 1])
        data = data[order]
        for mi, mol in enumerate(mol_ids):
            sel = data[data[:, 1] == mol]
            if len(sel) == 0:
                raise ValueError(f"Molecule id {mol} missing from frame {fi}.")
            m = np.array([mass_map.get(int(t), None) for t in sel[:, 2]])
            if any(x is None for x in m):
                bad = sel[[x is None for x in m], 2]
                raise ValueError(f"Unknown atom type(s) {bad} not in MASS_MAP -- edit SETTINGS.")
            com[fi, mi, :] = (sel[:, 3:6] * m[:, None]).sum(axis=0) / m.sum()
    times_ps = (np.array(timesteps, dtype=float) - timesteps[0]) * timestep_fs / 1000.0
    return times_ps, mol_ids, com


def autocorr_fft(x):
    """Autocorrelation via FFT (used inside the fast MSD algorithm)."""
    N = len(x)
    F = np.fft.fft(x, n=2 * N)
    PSD = F * F.conjugate()
    res = np.fft.ifft(PSD)
    res = res[:N].real
    n = N - np.arange(N)
    return res / n


def msd_fft_1d(r):
    """
    MSD(m) = <(r(k+m) - r(k))^2> averaged over all valid time origins k,
    computed in O(N log N) instead of the naive O(N^2) double loop.
    """
    N = len(r)
    D = np.square(r)
    D = np.append(D, 0)
    S2 = autocorr_fft(r)
    Q = 2 * D.sum()
    S1 = np.zeros(N)
    for m in range(N):
        Q -= D[m - 1] + D[N - m]
        S1[m] = Q / (N - m)
    return S1 - 2 * S2


def msd_fft_3d(traj):
    """traj: (n_frames, 3) array for one molecule -> per-axis + total MSD."""
    msd_x = msd_fft_1d(traj[:, 0])
    msd_y = msd_fft_1d(traj[:, 1])
    msd_z = msd_fft_1d(traj[:, 2])
    return msd_x, msd_y, msd_z, msd_x + msd_y + msd_z


def local_loglog_slope(t, msd, window=5):
    """Local slope of log(MSD) vs log(t) via a sliding linear fit."""
    slopes = np.full(len(t), np.nan)
    logt = np.log(t[1:])
    logmsd = np.log(np.clip(msd[1:], 1e-12, None))
    for i in range(window, len(logt) - window):
        xw = logt[i - window:i + window]
        yw = logmsd[i - window:i + window]
        A = np.vstack([xw, np.ones_like(xw)]).T
        s, _ = np.linalg.lstsq(A, yw, rcond=None)[0]
        slopes[i + 1] = s
    return slopes


def find_diffusive_window(t, msd, slope_tol, skip_fraction):
    """Find the longest contiguous stretch where local log-log slope ~ 1."""
    slopes = local_loglog_slope(t, msd)
    start_idx = int(len(t) * skip_fraction)
    ok = np.abs(slopes - 1.0) < slope_tol
    ok[:start_idx] = False

    best_start, best_len = None, 0
    cur_start, cur_len = None, 0
    for i, flag in enumerate(ok):
        if flag:
            if cur_start is None:
                cur_start = i
            cur_len += 1
            if cur_len > best_len:
                best_len, best_start = cur_len, cur_start
        else:
            cur_start, cur_len = None, 0

    if best_start is None or best_len < 10:
        # Fallback: no clean diffusive window found by slope criterion --
        # use the latter half of the (post-skip) trajectory instead, and
        # flag this clearly so it gets reported rather than silently used.
        lo = max(start_idx, len(t) // 2)
        return lo, len(t) - 1, False, slopes
    return best_start, best_start + best_len, True, slopes


def linear_fit(t, msd):
    """Unweighted least-squares fit of msd = a + b*t. Returns slope, intercept."""
    A = np.vstack([t, np.ones_like(t)]).T
    (b, a), *_ = np.linalg.lstsq(A, msd, rcond=None)
    return b, a


def analyze_segment(times_ps, com_segment):
    """
    Full per-segment analysis: per-molecule MSD, ensemble average,
    diffusive window detection, and D from isotropic + per-axis fits.
    com_segment: (n_frames, n_mol, 3), times_ps: (n_frames,) starting at 0.
    Returns a dict of results, or None if the segment is too short.
    """
    n_frames, n_mol, _ = com_segment.shape
    if n_frames < 30:
        return None

    msd_x_all = np.zeros((n_mol, n_frames))
    msd_y_all = np.zeros((n_mol, n_frames))
    msd_z_all = np.zeros((n_mol, n_frames))
    msd_tot_all = np.zeros((n_mol, n_frames))

    for mi in range(n_mol):
        mx, my, mz, mt = msd_fft_3d(com_segment[:, mi, :])
        msd_x_all[mi], msd_y_all[mi], msd_z_all[mi], msd_tot_all[mi] = mx, my, mz, mt

    msd_x = msd_x_all.mean(axis=0)
    msd_y = msd_y_all.mean(axis=0)
    msd_z = msd_z_all.mean(axis=0)
    msd_tot = msd_tot_all.mean(axis=0)
    msd_tot_std = msd_tot_all.std(axis=0)

    lo, hi, found_by_slope, slopes = find_diffusive_window(
        times_ps, msd_tot, SLOPE_TOL, SKIP_FRACTION
    )

    t_fit = times_ps[lo:hi]
    b_tot, a_tot = linear_fit(t_fit, msd_tot[lo:hi])
    b_x, _ = linear_fit(t_fit, msd_x[lo:hi])
    b_y, _ = linear_fit(t_fit, msd_y[lo:hi])
    b_z, _ = linear_fit(t_fit, msd_z[lo:hi])

    # D = slope / (2 * dimensionality); 3D isotropic total uses slope/6
    A2_PER_PS_TO_CM2_PER_S = 1e-4
    D_iso = (b_tot / 6.0) * A2_PER_PS_TO_CM2_PER_S
    D_x = (b_x / 2.0) * A2_PER_PS_TO_CM2_PER_S
    D_y = (b_y / 2.0) * A2_PER_PS_TO_CM2_PER_S
    D_z = (b_z / 2.0) * A2_PER_PS_TO_CM2_PER_S

    return dict(
        times_ps=times_ps, msd_x=msd_x, msd_y=msd_y, msd_z=msd_z,
        msd_tot=msd_tot, msd_tot_std=msd_tot_std, msd_tot_all=msd_tot_all,
        slopes=slopes, fit_lo=lo, fit_hi=hi, found_by_slope=found_by_slope,
        D_iso=D_iso, D_x=D_x, D_y=D_y, D_z=D_z,
        fit_t_start=times_ps[lo], fit_t_end=times_ps[hi - 1],
    )


def main():
    out = Path(OUT_DIR)
    out.mkdir(exist_ok=True)

    print(f"Parsing {DUMP_FILE} ...")
    timesteps, frames = parse_dump(DUMP_FILE)
    n_frames = len(frames)
    print(f"  {n_frames} frames, {frames[0].shape[0]} atoms per frame")

    times_ps, mol_ids, com = compute_com_trajectories(
        timesteps, frames, MASS_MAP, TIMESTEP_FS
    )
    n_mol = len(mol_ids)
    total_time_ps = times_ps[-1]
    print(f"  {n_mol} CO2 molecules, {total_time_ps:.1f} ps total trajectory")

    # ---- Full-trajectory analysis --------------------------------------
    full = analyze_segment(times_ps, com)
    if full is None:
        raise RuntimeError("Trajectory too short for meaningful MSD analysis.")

    # ---- Block-averaged uncertainty ------------------------------------
    block_len = n_frames // N_BLOCKS
    block_results = []
    for b in range(N_BLOCKS):
        s, e = b * block_len, (b + 1) * block_len
        if e - s < 30:
            continue
        seg_com = com[s:e]
        seg_t = times_ps[s:e] - times_ps[s]
        res = analyze_segment(seg_t, seg_com)
        if res is not None:
            block_results.append(res)

    D_iso_blocks = np.array([r["D_iso"] for r in block_results])
    D_x_blocks = np.array([r["D_x"] for r in block_results])
    D_y_blocks = np.array([r["D_y"] for r in block_results])
    D_z_blocks = np.array([r["D_z"] for r in block_results])

    # ---- Per-molecule trapping check ------------------------------------
    # Flag molecules whose individual total MSD at the end of the fit
    # window is far below the ensemble mean (possible cage-trapped guest).
    end_idx = full["fit_hi"] - 1
    per_mol_end_msd = full["msd_tot_all"][:, end_idx]
    ens_mean = per_mol_end_msd.mean()
    ens_std = per_mol_end_msd.std()
    trapped_flag = per_mol_end_msd < (ens_mean - 1.5 * ens_std)
    trapped_ids = mol_ids[trapped_flag]

    # =====================================================================
    # Plots
    # =====================================================================

    # 1. MSD vs t (linear) with fit window highlighted
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(full["times_ps"], full["msd_tot"], label="Ensemble-avg total MSD")
    ax.fill_between(
        full["times_ps"],
        full["msd_tot"] - full["msd_tot_std"],
        full["msd_tot"] + full["msd_tot_std"],
        alpha=0.2, label="+/- 1 std across molecules"
    )
    t_fit = full["times_ps"][full["fit_lo"]:full["fit_hi"]]
    b_tot = full["D_iso"] * 6.0 / 1e-4
    a_tot = full["msd_tot"][full["fit_lo"]] - b_tot * t_fit[0]
    ax.plot(t_fit, a_tot + b_tot * t_fit, "r--", lw=2, label="Linear fit (diffusive window)")
    ax.axvspan(full["fit_t_start"], full["fit_t_end"], color="red", alpha=0.08)
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel(r"MSD ($\mathrm{\AA}^2$)")
    ax.set_title("CO2 center-of-mass MSD in MOF-5")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "msd_vs_time.png", dpi=150)
    plt.close(fig)

    # 2. log-log MSD with local slope, to justify the fit window
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 8), sharex=True)
    ax1.loglog(full["times_ps"][1:], full["msd_tot"][1:])
    ax1.axvspan(full["fit_t_start"], full["fit_t_end"], color="red", alpha=0.08)
    ax1.set_ylabel(r"MSD ($\mathrm{\AA}^2$)")
    ax1.set_title("Diffusive-regime identification (log-log slope -> 1 = Fickian)")
    ax2.plot(full["times_ps"], full["slopes"])
    ax2.axhline(1.0, color="k", lw=1, ls=":")
    ax2.axhspan(1 - SLOPE_TOL, 1 + SLOPE_TOL, color="green", alpha=0.1)
    ax2.axvspan(full["fit_t_start"], full["fit_t_end"], color="red", alpha=0.08)
    ax2.set_xlabel("Time (ps)")
    ax2.set_ylabel("Local d(log MSD)/d(log t)")
    ax2.set_xscale("log")
    fig.tight_layout()
    fig.savefig(out / "diffusive_regime_check.png", dpi=150)
    plt.close(fig)

    # 3. Per-axis MSD (anisotropy check)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(full["times_ps"], full["msd_x"], label="x")
    ax.plot(full["times_ps"], full["msd_y"], label="y")
    ax.plot(full["times_ps"], full["msd_z"], label="z")
    ax.axvspan(full["fit_t_start"], full["fit_t_end"], color="red", alpha=0.08)
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel(r"MSD component ($\mathrm{\AA}^2$)")
    ax.set_title("Per-axis MSD (anisotropy check, cubic MOF-5 expects x~y~z)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "msd_per_axis.png", dpi=150)
    plt.close(fig)

    # 4. All individual molecule MSD curves (trapping / outlier check)
    fig, ax = plt.subplots(figsize=(7, 5))
    for mi in range(n_mol):
        color = "tab:red" if trapped_flag[mi] else "tab:gray"
        alpha = 0.9 if trapped_flag[mi] else 0.35
        ax.plot(full["times_ps"], full["msd_tot_all"][mi], color=color, alpha=alpha, lw=1)
    ax.plot(full["times_ps"], full["msd_tot"], color="k", lw=2, label="Ensemble average")
    ax.axvspan(full["fit_t_start"], full["fit_t_end"], color="red", alpha=0.08)
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel(r"MSD ($\mathrm{\AA}^2$)")
    ax.set_title("Per-molecule MSD (red = possible cage-trapped outlier)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "msd_per_molecule.png", dpi=150)
    plt.close(fig)

    # =====================================================================
    # Data export
    # =====================================================================
    np.savetxt(
        out / "msd_data.csv",
        np.column_stack([full["times_ps"], full["msd_tot"], full["msd_tot_std"],
                          full["msd_x"], full["msd_y"], full["msd_z"]]),
        header="time_ps,msd_total_A2,msd_total_std_A2,msd_x_A2,msd_y_A2,msd_z_A2",
        delimiter=",", comments=""
    )

    # =====================================================================
    # Summary
    # =====================================================================
    lines = []
    lines.append("=" * 70)
    lines.append("CO2 self-diffusion in MOF-5 -- analysis summary")
    lines.append("=" * 70)
    lines.append("")
    lines.append("SYSTEM")
    lines.append(f"  CO2 molecules analyzed          : {n_mol}")
    lines.append(f"  Framework                       : UFF, fully frozen (no framework dynamics)")
    lines.append(f"  CO2 model                       : TraPPE, rigid body (fix rigid/nve/small)")
    lines.append(f"  Production ensemble              : NVE (no thermostat)")
    lines.append(f"  Trajectory length                : {total_time_ps:.1f} ps ({n_frames} frames, "
                 f"{DUMP_EVERY} steps / {TIMESTEP_FS} fs per frame)")
    lines.append("")
    lines.append("DIFFUSIVE WINDOW")
    if full["found_by_slope"]:
        lines.append(f"  Identified automatically from local log-log slope ~ 1 +/- {SLOPE_TOL}")
    else:
        lines.append("  WARNING: no window satisfied the slope~1 criterion automatically;")
        lines.append("  fell back to the latter half of the (post-transient) trajectory.")
        lines.append("  Inspect diffusive_regime_check.png before trusting the reported D --")
        lines.append("  you likely need a longer production run.")
    lines.append(f"  Fit window                       : {full['fit_t_start']:.1f} - {full['fit_t_end']:.1f} ps")
    lines.append("")
    lines.append("SELF-DIFFUSION COEFFICIENT (full-trajectory fit)")
    lines.append(f"  D_isotropic  (3D, total MSD/6t)  : {full['D_iso']:.3e} cm^2/s")
    lines.append(f"  D_x (MSD_x/2t)                   : {full['D_x']:.3e} cm^2/s")
    lines.append(f"  D_y (MSD_y/2t)                    : {full['D_y']:.3e} cm^2/s")
    lines.append(f"  D_z (MSD_z/2t)                    : {full['D_z']:.3e} cm^2/s")
    lines.append("")
    lines.append("BLOCK-AVERAGED UNCERTAINTY")
    lines.append(f"  Number of blocks used             : {len(block_results)} / {N_BLOCKS} requested")
    if len(D_iso_blocks) > 1:
        lines.append(f"  D_isotropic (mean +/- std)         : "
                      f"{D_iso_blocks.mean():.3e} +/- {D_iso_blocks.std():.3e} cm^2/s")
        lines.append(f"  D_x (mean +/- std)                 : "
                      f"{D_x_blocks.mean():.3e} +/- {D_x_blocks.std():.3e} cm^2/s")
        lines.append(f"  D_y (mean +/- std)                 : "
                      f"{D_y_blocks.mean():.3e} +/- {D_y_blocks.std():.3e} cm^2/s")
        lines.append(f"  D_z (mean +/- std)                 : "
                      f"{D_z_blocks.mean():.3e} +/- {D_z_blocks.std():.3e} cm^2/s")
        lines.append("  NOTE: block std here reflects statistical scatter across time")
        lines.append("  windows of ONE trajectory, not independent replicate runs. For")
        lines.append("  full Q1-level rigor, running 2-3 replicate production runs from")
        lines.append("  different initial velocity seeds and combining D across replicates")
        lines.append("  AND blocks is the more defensible uncertainty estimate.")
    else:
        lines.append("  Not enough blocks produced a valid diffusive window -- trajectory")
        lines.append("  may be too short relative to N_BLOCKS. Consider a longer run or")
        lines.append("  fewer blocks.")
    lines.append("")
    lines.append("ANISOTROPY CHECK")
    dvals = np.array([full["D_x"], full["D_y"], full["D_z"]])
    spread = (dvals.max() - dvals.min()) / dvals.mean() * 100
    lines.append(f"  Spread among Dx, Dy, Dz            : {spread:.1f}% of mean")
    if spread > 30:
        lines.append("  WARNING: MOF-5 has cubic symmetry, so a well-converged, adequately")
        lines.append("  sampled diffusion trajectory should show Dx ~ Dy ~ Dz. This large a")
        lines.append("  spread suggests insufficient sampling (too few molecules/too short a")
        lines.append("  run for the per-axis statistics to converge), not necessarily real")
        lines.append("  physical anisotropy. Longer/replicate runs are recommended before")
        lines.append("  reporting the isotropic D with confidence.")
    else:
        lines.append("  Spread is modest and consistent with expected cubic-symmetry isotropy.")
    lines.append("")
    lines.append("PER-MOLECULE TRAPPING CHECK")
    if len(trapped_ids) > 0:
        lines.append(f"  {len(trapped_ids)} of {n_mol} molecules show MSD far below the ensemble")
        lines.append(f"  mean at the end of the fit window (molecule IDs: {list(trapped_ids)}).")
        lines.append("  This can indicate cage-trapped/low-mobility guests contributing little")
        lines.append("  to net diffusion -- inspect msd_per_molecule.png. Consider whether to")
        lines.append("  report a trimmed-mean D excluding outliers, or discuss this")
        lines.append("  heterogeneity explicitly in the paper (it may be physically real).")
    else:
        lines.append("  No molecules flagged as outliers/trapped relative to the ensemble.")
    lines.append("")
    lines.append("FILES WRITTEN")
    lines.append("  msd_vs_time.png            - MSD(t) with fit window and linear fit overlaid")
    lines.append("  diffusive_regime_check.png - log-log MSD + local slope, justifies fit window")
    lines.append("  msd_per_axis.png           - per-axis MSD, anisotropy check")
    lines.append("  msd_per_molecule.png       - all 20 individual-molecule MSD curves")
    lines.append("  msd_data.csv               - raw ensemble MSD(t) data (time, total, x, y, z)")
    lines.append("  summary.txt                - this file")
    lines.append("")
    lines.append("CAVEATS FOR THE METHODS SECTION")
    lines.append("  - Framework was treated as fully rigid/frozen; framework flexibility")
    lines.append("    (breathing modes) can measurably affect diffusion coefficients in MOFs")
    lines.append("    and should be mentioned as a modeling choice/limitation if not also")
    lines.append("    tested with a flexible framework for comparison.")
    lines.append("  - This is a single production trajectory. Reviewers will likely expect")
    lines.append("    either multiple independent replicates (different initial velocities)")
    lines.append("    or a substantially longer single run, with uncertainty reported across")
    lines.append("    replicates rather than only across time-blocks of one run.")
    lines.append("  - D reported here is uncorrected for finite-size effects; note this if")
    lines.append("    comparing directly to experimental or other simulation literature values.")

    summary_text = "\n".join(lines)
    (out / "summary.txt").write_text(summary_text)
    print(summary_text)
    print(f"\nAll outputs written to: {out.resolve()}")


if __name__ == "__main__":
    main()
