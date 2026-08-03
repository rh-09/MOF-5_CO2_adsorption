#!/usr/bin/env python3
"""
extract_equilibrium_frame.py

Purpose
-------
1. Read a RASPA2 .data output file, confirm the run has reached equilibrium
   (by comparing the cumulative average loading at the run's midpoint vs.
   its final value), and determine the target integer equilibrium loading.
2. Read the corresponding per-component movie .pdb file, find the frame(s)
   whose molecule count is closest to that target loading (ties broken by
   taking the LATEST such frame), and extract that single frame to its own
   .pdb file.

Usage
-----
Edit the three filenames in the CONFIG section below, then run:

    python3 extract_equilibrium_frame.py

Notes
-----
- Assumes CO2 (3 atoms/molecule) in the movie file. Change ATOMS_PER_MOLECULE
  if you're working with a different adsorbate.
- Equilibration tolerance is a standard block-averaging check: if the
  cumulative average at the ~50% mark agrees with the final cumulative
  average within EQUILIBRATION_TOLERANCE (relative), the run is considered
  equilibrated. This is a common, simple convergence diagnostic for GCMC.
"""

import re
import sys

# ---------------------------------------------------------------------------
# CONFIG -- edit these for your files
# ---------------------------------------------------------------------------
DATA_FILE = r"C:\Users\User\Desktop\ex\Output\System_0\output_MOF-5_raspa_1.1.1_298.000000_100000.data"   # RASPA2 data output file
MOVIE_FILE = r"C:\Users\User\Desktop\ex\Movies\System_0\Movie_MOF-5_raspa_1.1.1_298.000000_100000.000000_component_CO2_0.pdb"   # RASPA2 per-component movie file
OUTPUT_FRAME_FILE = "best_frame_CO2.pdb"   # where the extracted frame is written

ATOMS_PER_MOLECULE = 3          # CO2 = 3 atoms (1 C + 2 O)
EQUILIBRATION_TOLERANCE = 0.02  # 2% relative difference tolerance
COMPONENT_NAME_HINT = "CO2"     # used only in printed messages


# ---------------------------------------------------------------------------
# Step 1: parse the .data output file
# ---------------------------------------------------------------------------
def parse_data_file(path):
    """
    Extract (cycle, cumulative_average_loading) pairs from a RASPA2 .data
    output file by matching:
        Current cycle: <N> out of <M>
        ... current number of integer/fractional/reaction molecules: X/Y/Z (avg. AVG)
    """
    cycle_re = re.compile(r"Current cycle:\s*(\d+)\s*out of\s*(\d+)")
    avg_re = re.compile(
        r"current number of integer/fractional/reaction molecules:\s*"
        r"\d+/\d+/\d+\s*\(avg\.\s*([\d.]+)\)"
    )

    entries = []  # list of (cycle, avg_loading)
    current_cycle = None

    with open(path, "r") as f:
        for line in f:
            cyc_match = cycle_re.search(line)
            if cyc_match:
                current_cycle = int(cyc_match.group(1))
                continue

            avg_match = avg_re.search(line)
            if avg_match and current_cycle is not None:
                avg_loading = float(avg_match.group(1))
                entries.append((current_cycle, avg_loading))
                current_cycle = None  # reset until next "Current cycle" line

    if not entries:
        sys.exit(
            f"ERROR: No loading data found in '{path}'. "
            "Check the file path and format."
        )

    return entries


def check_equilibration(entries):
    """
    Standard block-averaging equilibration check: compare the cumulative
    average loading at the run's midpoint cycle against the final cumulative
    average. If they agree within EQUILIBRATION_TOLERANCE (relative), the
    run is considered equilibrated.
    """
    final_cycle, final_avg = entries[-1]
    mid_cycle_target = final_cycle / 2.0

    # find the entry whose cycle is closest to the midpoint
    mid_entry = min(entries, key=lambda e: abs(e[0] - mid_cycle_target))
    mid_cycle, mid_avg = mid_entry

    rel_diff = abs(final_avg - mid_avg) / final_avg

    print("=" * 60)
    print("EQUILIBRATION CHECK")
    print("=" * 60)
    print(f"Midpoint  -> cycle {mid_cycle:>7d}, cumulative avg = {mid_avg:.5f}")
    print(f"Final     -> cycle {final_cycle:>7d}, cumulative avg = {final_avg:.5f}")
    print(f"Relative difference: {rel_diff * 100:.3f}%  "
          f"(tolerance: {EQUILIBRATION_TOLERANCE * 100:.1f}%)")

    equilibrated = rel_diff <= EQUILIBRATION_TOLERANCE

    if equilibrated:
        print("Result: EQUILIBRATED (plateau confirmed).")
    else:
        print("Result: NOT clearly equilibrated -- the average is still "
              "drifting. Consider extending the run before trusting this "
              "target loading.")
    print("=" * 60)

    return equilibrated, final_avg


# ---------------------------------------------------------------------------
# Step 2: parse the movie file and find the closest frame
# ---------------------------------------------------------------------------
def parse_movie_file(path, atoms_per_molecule):
    """
    Parse a RASPA2 per-component movie .pdb file (MODEL/ATOM/ENDMDL blocks)
    and return a list of (model_number, molecule_count, start_line, end_line)
    for every frame, plus the full list of raw lines for later extraction.
    """
    with open(path, "r") as f:
        lines = f.readlines()

    frames = []  # (model_number, molecule_count, start_idx, end_idx)
    current_model = None
    atom_count = 0
    start_idx = None

    for idx, line in enumerate(lines):
        if line.startswith("MODEL"):
            current_model = int(line.split()[1])
            atom_count = 0
            start_idx = idx
        elif line.startswith("ATOM"):
            atom_count += 1
        elif line.startswith("ENDMDL"):
            if current_model is not None:
                molecule_count = atom_count / atoms_per_molecule
                frames.append((current_model, molecule_count, start_idx, idx))
            current_model = None

    if not frames:
        sys.exit(f"ERROR: No MODEL/ENDMDL frames found in '{path}'.")

    return frames, lines


def find_closest_frame(frames, target):
    """
    Find the frame(s) whose molecule count is closest to `target`.
    Ties are broken by taking the LATEST frame (highest model number).
    """
    min_diff = min(abs(mc - target) for (_, mc, _, _) in frames)
    candidates = [fr for fr in frames if abs(fr[1] - target) == min_diff]
    best = max(candidates, key=lambda fr: fr[0])  # latest model number wins ties
    return best, candidates


def extract_frame(lines, start_idx, end_idx, output_path):
    """Write lines[start_idx:end_idx+1] (inclusive of MODEL/ENDMDL) to a file."""
    with open(output_path, "w") as f:
        f.writelines(lines[start_idx:end_idx + 1])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print(f"\nReading RASPA2 output data file: {DATA_FILE}")
    entries = parse_data_file(DATA_FILE)

    equilibrated, final_avg = check_equilibration(entries)

    target_loading = round(final_avg)
    print(f"\nTarget equilibrium loading (rounded to nearest integer): "
          f"{target_loading} molecules of {COMPONENT_NAME_HINT}")

    if not equilibrated:
        print("\nWARNING: proceeding anyway, but treat the extracted frame "
              "with caution since equilibration was not clearly confirmed.")

    print(f"\nReading movie file: {MOVIE_FILE}")
    frames, lines = parse_movie_file(MOVIE_FILE, ATOMS_PER_MOLECULE)
    print(f"Found {len(frames)} frames in movie file.")

    best_frame, tied_candidates = find_closest_frame(frames, target_loading)
    best_model, best_count, best_start, best_end = best_frame

    print("\n" + "=" * 60)
    print("FRAME SELECTION")
    print("=" * 60)
    if len(tied_candidates) > 1:
        tied_models = [c[0] for c in tied_candidates]
        print(f"{len(tied_candidates)} frames tied at the closest distance "
              f"(models: {tied_models}). Selecting the latest: model {best_model}.")
    print(f"Selected frame -> MODEL {best_model}, "
          f"{best_count:.0f} molecules "
          f"(target was {target_loading})")
    print("=" * 60)

    extract_frame(lines, best_start, best_end, OUTPUT_FRAME_FILE)
    print(f"\nExtracted frame written to: {OUTPUT_FRAME_FILE}")


if __name__ == "__main__":
    main()
