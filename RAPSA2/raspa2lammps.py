#!/usr/bin/env python3
"""
raspa2lammps.py — merge a single extracted RASPA2 movie frame (an
equilibrium-representative GCMC-loaded configuration, e.g. from
extract_equilibrium_frame.py) with a LAMMPS framework data file (from
lammps-interface) into one complete LAMMPS data file, for any MOF and any
adsorbate.

WORKFLOW THIS SUPPORTS
-----------------------
1. Build the empty-framework LAMMPS files with lammps-interface (data.<MOF>).
2. Run the matching RASPA2 GCMC simulation (same CIF, same UnitCells).
3. Confirm equilibration and extract the frame whose loading is closest to
   the equilibrium average (e.g. with extract_equilibrium_frame.py) — a
   single-MODEL PDB file containing only that adsorbate species' atoms.
4. Run this script pointing at:
     - the framework data file
     - the extracted frame PDB file
     - the RASPA working directory (pseudo_atoms.def, force_field_mixing_rules.def,
       and every <MoleculeName>.def used in the simulation)

USAGE
-----
    python3 raspa2lammps.py \
        --framework data.MOF-5 \
        --frame best_frame_CO2.pdb \
        --species CO2 \
        --raspa-dir /path/to/raspa/working/directory \
        --out data.MOF-5_loaded

If --species is omitted, the script tries to auto-detect it by matching a
.def filename in --raspa-dir against the frame file's own filename (e.g.
"best_frame_CO2.pdb" -> "CO2"). If that's ambiguous, pass --species explicitly.

WHAT IT DOES NOT DO
--------------------
- It does not run RASPA, lammps-interface, or the frame-extraction step for you.
- It assumes a RIGID framework (the overwhelmingly common case for GCMC
  adsorption screening). Extracted movie frames only ever contain the
  adsorbate species' atoms, never framework atoms, so there is no
  flexible-framework option here — framework coordinates always come from
  the framework data file itself.
- Bond connectivity for each adsorbate is taken from the "# Bond stretch"
  section of its .def file. If a species has no bonds defined there
  (e.g. a single-atom species like Ar/Kr/Xe, or a rigid molecule whose
  .def never declares bonds), no Bonds are written for it — you'll need
  `fix rigid/small molecule` (or similar) in LAMMPS for those.
- Assumes the extracted frame lists each molecule's atoms in the same
  order used throughout the molecule's own .def file (index 0, 1, 2, ...) —
  this is RASPA2's standard internal ordering and is consistent between
  its movie output and its .def-file atom declarations.
"""

import argparse
import glob
import math
import os
import re
import sys

KB_KCAL = 0.0019872041  # Boltzmann constant, kcal/(mol*K) - converts RASPA epsilon (K) to kcal/mol


# ============================================================================
# Parsers
# ============================================================================

def parse_pseudo_atoms(path):
    """Return {type_name: {'mass': float, 'charge': float}}"""
    out = {}
    with open(path, "r", newline=None) as f:
        lines = [l.rstrip("\n") for l in f]
    started = False
    for l in lines:
        s = l.strip()
        if not s or s.startswith("#"):
            continue
        if s.lower().startswith("number of pseudo atoms") or re.match(r"^\d+$", s):
            started = True
            continue
        parts = s.split()
        if len(parts) < 7:
            continue
        name = parts[0]
        try:
            mass = float(parts[5])
            charge = float(parts[6])
        except ValueError:
            continue
        out[name] = {"mass": mass, "charge": charge}
    return out


def parse_mixing_rules(path):
    """Return {type_name: {'style': str, 'eps_K': float or None, 'sigma': float or None}}"""
    out = {}
    with open(path, "r", newline=None) as f:
        for l in f:
            s = l.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) < 2:
                continue
            name, style = parts[0], parts[1].lower()
            if style == "lennard-jones" and len(parts) >= 4:
                out[name] = {"style": style, "eps_K": float(parts[2]), "sigma": float(parts[3])}
            elif style == "none":
                out[name] = {"style": style, "eps_K": None, "sigma": None}
            else:
                out[name] = {"style": style, "eps_K": None, "sigma": None}
    return out


def parse_molecule_def(path):
    """
    Parse a RASPA molecule .def file (CO2.def-style).
    Returns:
      n_atoms: total atom count declared at top of file
      atom_type_by_index: {global_atom_index: pseudo_atom_name}
      bonds: [(atom1_idx, atom2_idx, bond_keyword)]
    """
    with open(path, "r", newline=None) as f:
        raw_lines = [l.rstrip("\n").rstrip("\r") for l in f]
    # Keep only non-blank lines; comments are NOT stripped here because the
    # RASPA parser itself treats "the next line" positionally, not by content -
    # we replicate that same positional logic to stay faithful to the format.
    lines = [l for l in raw_lines if l.strip() != ""]

    idx = 0
    def next_line():
        nonlocal idx
        l = lines[idx]
        idx += 1
        return l.strip()

    next_line()  # comment: critical constants header
    next_line()  # Tc
    next_line()  # Pc
    next_line()  # acentric factor
    next_line()  # comment: number of atoms
    n_atoms = int(next_line().split()[0])
    next_line()  # comment: number of groups
    n_groups = int(next_line().split()[0])

    atom_type_by_index = {}
    for _ in range(n_groups):
        next_line()  # comment (e.g. "# rigid group")
        next_line()  # keyword: rigid / flexible
        next_line()  # comment: number of atoms
        group_atom_count = int(next_line().split()[0])
        if group_atom_count > 0:
            next_line()  # comment: atomic positions
            for _ in range(group_atom_count):
                toks = next_line().split()
                gidx, name = int(toks[0]), toks[1]
                atom_type_by_index[gidx] = name

    next_line()  # comment: Chiral centers Bond BondDipoles ...
    counts = next_line().split()
    n_bonds_declared = int(counts[1]) if len(counts) > 1 else 0

    bonds = []
    if n_bonds_declared > 0:
        next_line()  # comment: "# Bond stretch: ..."
        for _ in range(n_bonds_declared):
            toks = next_line().split()
            a1, a2, kw = int(toks[0]), int(toks[1]), toks[2]
            bonds.append((a1, a2, kw))

    return {"n_atoms": n_atoms, "atom_type_by_index": atom_type_by_index, "bonds": bonds}


def find_molecule_def_files(raspa_dir):
    """
    Auto-discover <Name>.def molecule files in the RASPA working directory,
    excluding the reserved force-field filenames.
    Returns {ComponentName: path}
    """
    reserved = {"pseudo_atoms.def", "force_field.def", "force_field_mixing_rules.def"}
    out = {}
    for p in glob.glob(os.path.join(raspa_dir, "*.def")):
        base = os.path.basename(p)
        if base in reserved:
            continue
        name = os.path.splitext(base)[0]
        out[name] = p
    return out


def parse_frame_pdb(path):
    """
    Parse a single extracted movie-frame PDB (one MODEL/ENDMDL block
    containing only one adsorbate species' atoms, as produced by
    extract_equilibrium_frame.py).

    Returns:
      cell_vectors: 3x3 list (a, b, c vectors), derived from the CRYST1
                    record using the standard crystallographic->Cartesian
                    convention (a along x, b in the xy-plane) -- the same
                    convention RASPA2/LAMMPS triclinic boxes use. None if
                    no CRYST1 record is present.
      atom_coords:  ordered list of (x, y, z) for every ATOM line, in file
                    order (this order is what lets us regroup atoms into
                    molecules below).
    """
    cell_vectors = None
    atom_coords = []
    with open(path, "r", newline=None) as f:
        for line in f:
            if line.startswith("CRYST1"):
                parts = line.split()
                a, b, c = float(parts[1]), float(parts[2]), float(parts[3])
                alpha, beta, gamma = float(parts[4]), float(parts[5]), float(parts[6])
                al, be, ga = math.radians(alpha), math.radians(beta), math.radians(gamma)
                ax, ay, az = a, 0.0, 0.0
                bx, by, bz = b * math.cos(ga), b * math.sin(ga), 0.0
                cx = c * math.cos(be)
                cy = c * (math.cos(al) - math.cos(be) * math.cos(ga)) / math.sin(ga)
                cz_sq = c * c - cx * cx - cy * cy
                cz = math.sqrt(cz_sq) if cz_sq > 0 else 0.0
                cell_vectors = [(ax, ay, az), (bx, by, bz), (cx, cy, cz)]
            elif line.startswith("ATOM"):
                toks = line.split()
                x, y, z = float(toks[4]), float(toks[5]), float(toks[6])
                atom_coords.append((x, y, z))

    if not atom_coords:
        sys.exit(f"ERROR: no ATOM lines found in frame file '{path}'.")

    return cell_vectors, atom_coords


def build_restart_like_from_frame(frame_path, species_name, n_atoms_per_molecule):
    """
    Regroup a single extracted frame's flat atom list into per-molecule
    positions, and package it into the same shape the merge logic below
    expects from a restart file:
      cell_vectors: 3x3 list (or None)
      components: {0: {'name': species_name, 'n_molecules': int}}
      adsorbate_positions: {(0, mol, atom): (x, y, z)}
      cation_positions: {}   (frame files never contain cations)
      framework_positions: {} (frame files never contain framework atoms --
                               framework coordinates always come from the
                               framework data file itself)
    """
    cell_vectors, atom_coords = parse_frame_pdb(frame_path)

    if len(atom_coords) % n_atoms_per_molecule != 0:
        sys.exit(
            f"ERROR: frame file '{frame_path}' has {len(atom_coords)} atoms, "
            f"which is not evenly divisible by {n_atoms_per_molecule} "
            f"atoms/molecule declared for species '{species_name}'. "
            f"Check --species and the molecule's .def file."
        )

    n_mol = len(atom_coords) // n_atoms_per_molecule
    adsorbate_positions = {}
    for mol in range(n_mol):
        for local_atom in range(n_atoms_per_molecule):
            xyz = atom_coords[mol * n_atoms_per_molecule + local_atom]
            adsorbate_positions[(0, mol, local_atom)] = xyz

    components = {0: {"name": species_name, "n_molecules": n_mol}}

    return {
        "cell_vectors": cell_vectors,
        "components": components,
        "adsorbate_positions": adsorbate_positions,
        "cation_positions": {},
        "framework_positions": {},
    }


def autodetect_species_from_filename(frame_path, mol_def_files):
    """
    Best-effort guess at the adsorbate species name by checking whether any
    known .def filename (from --raspa-dir) appears in the frame file's own
    filename (e.g. 'best_frame_CO2.pdb' or
    'Movie_..._component_CO2_0.pdb' -> 'CO2'). Returns None if zero or
    multiple .def names match, so the caller can fall back to requiring
    --species explicitly.
    """
    base = os.path.basename(frame_path)
    matches = [name for name in mol_def_files if name in base]
    if len(matches) == 1:
        return matches[0]
    return None


def parse_framework_lammps(path):
    with open(path, "r", newline=None) as f:
        raw = f.read()
    lines = raw.splitlines()

    header_counts = {}
    for key in ["atoms", "bonds", "angles", "dihedrals", "impropers",
                "atom types", "bond types", "angle types", "dihedral types", "improper types"]:
        for l in lines:
            m = re.match(rf"^\s*(\d+)\s+{key}\s*$", l)
            if m:
                header_counts[key] = int(m.group(1))
                break

    box_lines = [l for l in lines if re.search(r"xlo xhi|ylo yhi|zlo zhi|xy xz yz", l)]

    section_names = ["Masses", "Bond Coeffs", "Angle Coeffs", "Dihedral Coeffs",
                      "Improper Coeffs", "Pair Coeffs", "Atoms", "Bonds",
                      "Angles", "Dihedrals", "Impropers", "Velocities"]
    found = sorted((i, l.strip()) for i, l in enumerate(lines) if l.strip() in section_names)

    def get_block(name):
        for pos, (i, nm) in enumerate(found):
            if nm == name:
                start = i + 2
                end = found[pos + 1][0] if pos + 1 < len(found) else len(lines)
                return [x for x in lines[start:end] if x.strip()]
        return []

    blocks = {name: get_block(name) for name in section_names}

    return {"header_counts": header_counts, "box_lines": box_lines, "blocks": blocks}


# ============================================================================
# Box consistency check
# ============================================================================

def box_vectors_from_lammps(box_lines):
    xlo = xhi = ylo = yhi = zlo = zhi = xy = xz = yz = 0.0
    for l in box_lines:
        parts = l.split()
        if "xlo" in l:
            xlo, xhi = float(parts[0]), float(parts[1])
        elif "ylo" in l:
            ylo, yhi = float(parts[0]), float(parts[1])
        elif "zlo" in l:
            zlo, zhi = float(parts[0]), float(parts[1])
        elif "xy" in l:
            xy, xz, yz = float(parts[0]), float(parts[1]), float(parts[2])
    a = (xhi - xlo, 0.0, 0.0)
    b = (xy, yhi - ylo, 0.0)
    c = (xz, yz, zhi - zlo)
    return [a, b, c]


def check_box_consistency(frame_cell_vectors, lammps_box_lines, tol=0.05):
    lammps_vecs = box_vectors_from_lammps(lammps_box_lines)
    if frame_cell_vectors is None or None in frame_cell_vectors:
        print("WARNING: could not read cell vectors from frame file (no CRYST1 record?); skipping box consistency check.")
        return
    for label, rv, lv in zip("abc", frame_cell_vectors, lammps_vecs):
        diff = max(abs(rv[i] - lv[i]) for i in range(3))
        if diff > tol:
            print(f"WARNING: cell vector {label} differs between frame file and framework data file "
                  f"by {diff:.4f} A (frame={rv}, lammps={lv}).")
            print("         This usually means the framework data file was built with a different "
                  "UnitCells replication than the RASPA run. Coordinates will NOT line up correctly.")
    print("Box consistency check: frame vs framework box vectors compared (see warnings above, if any).")


def looks_like_lammps_data_file(path):
    """Heuristic: a LAMMPS data file has an 'N atoms' line within its first ~40 lines."""
    try:
        with open(path, "r", newline=None, errors="ignore") as f:
            for i, l in enumerate(f):
                if i > 40:
                    break
                if re.match(r"^\s*\d+\s+atoms\s*$", l):
                    return True
    except (IsADirectoryError, PermissionError):
        return False
    return False


def looks_like_frame_pdb_file(path):
    """Heuristic: an extracted single-frame PDB has both a MODEL and an ATOM line."""
    if path.endswith((".py", ".def", ".cif", ".input", ".pyc")):
        return False
    try:
        with open(path, "r", newline=None, errors="ignore") as f:
            head = f.read(4000)
    except (IsADirectoryError, PermissionError, UnicodeDecodeError):
        return False
    return ("MODEL" in head) and ("ATOM" in head) and ("cell-vector-a:" not in head)


def autodetect_framework(directory):
    """Find the LAMMPS framework data file in `directory`."""
    self_path = os.path.abspath(__file__)
    excluded_ext = (".py", ".def", ".cif", ".input", ".pyc", ".md", ".txt")
    candidates = []
    for p in sorted(glob.glob(os.path.join(directory, "*"))):
        if os.path.isdir(p) or os.path.abspath(p) == self_path:
            continue
        base = os.path.basename(p)
        if base.startswith(".") or base.endswith(excluded_ext):
            continue
        # Never pick up this script's own prior output (our naming convention
        # appends "_loaded"). Without this, re-running in the same folder
        # would merge CO2 into an already-loaded file a second time.
        if "_loaded" in base.lower():
            continue
        # fast filename hint first (LAMMPS convention: "data.<name>")
        if base.startswith("data.") or base.startswith("data_"):
            candidates.append(p)
        elif looks_like_lammps_data_file(p):
            candidates.append(p)
    candidates = sorted(set(candidates))
    if len(candidates) == 0:
        sys.exit(f"ERROR: could not auto-detect a LAMMPS framework data file in '{directory}'.\n"
                  f"       Expected a file like 'data.<MOF-name>'. Pass --framework explicitly if it's named differently.")
    if len(candidates) > 1:
        sys.exit(f"ERROR: multiple possible framework data files found and auto-detection can't safely "
                  f"choose between them:\n  " + "\n  ".join(candidates) +
                  f"\n       Re-run with --framework <path> to specify which one to use.")
    print(f"Auto-detected framework data file: {candidates[0]}")
    return candidates[0]


def autodetect_frame(directory):
    """Find an extracted single-frame adsorbate PDB file in `directory`."""
    self_path = os.path.abspath(__file__)
    excluded_ext = (".py", ".def", ".cif", ".input", ".pyc", ".md", ".txt")
    candidates = []
    for p in sorted(glob.glob(os.path.join(directory, "*"))):
        if os.path.isdir(p) or os.path.abspath(p) == self_path:
            continue
        base = os.path.basename(p)
        if base.endswith(excluded_ext):
            continue
        if "best_frame" in base.lower() or looks_like_frame_pdb_file(p):
            candidates.append(p)
    candidates = sorted(set(candidates))
    if len(candidates) == 0:
        sys.exit(f"ERROR: could not auto-detect an extracted frame PDB file in '{directory}'.\n"
                  f"       Expected something like 'best_frame_CO2.pdb' "
                  f"(from extract_equilibrium_frame.py). Pass --frame explicitly if needed.")
    if len(candidates) > 1:
        sys.exit(f"ERROR: multiple possible frame files found and auto-detection can't safely "
                  f"choose between them:\n  " + "\n  ".join(candidates) +
                  f"\n       Re-run with --frame <path> to specify which one to use.")
    print(f"Auto-detected frame file: {candidates[0]}")
    return candidates[0]


def autodetect_out_path(directory, framework_path):
    base = os.path.basename(framework_path)
    if base.startswith("data."):
        out_name = base + "_loaded"
    else:
        out_name = "data." + base + "_loaded"
    return os.path.join(directory, out_name)


# ============================================================================
# Main merge
# ============================================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default=".",
                     help="Directory to auto-detect everything from (default: current directory). "
                          "Use this if you don't pass --framework/--frame/--raspa-dir explicitly.")
    ap.add_argument("--framework", default=None, help="LAMMPS data file for the empty framework. Auto-detected from --dir if omitted.")
    ap.add_argument("--frame", default=None, help="Extracted single-frame adsorbate PDB file (e.g. best_frame_CO2.pdb, from extract_equilibrium_frame.py). Auto-detected from --dir if omitted.")
    ap.add_argument("--species", default=None, help="Adsorbate species name, must match a <species>.def file in --raspa-dir (e.g. CO2). Auto-detected from the --frame filename if omitted and unambiguous.")
    ap.add_argument("--raspa-dir", default=None, help="Directory containing pseudo_atoms.def, force_field_mixing_rules.def, and all <Molecule>.def files. Defaults to --dir.")
    ap.add_argument("--out", default=None, help="Output LAMMPS data file path. Auto-named from the framework file if omitted.")
    ap.add_argument("--bond-k", type=float, default=600.0,
                     help="Harmonic bond force constant (kcal/mol/A^2) used as a placeholder for RIGID_BOND entries (default: 600.0)")
    args = ap.parse_args()

    base_dir = os.path.abspath(args.dir)
    if not os.path.isdir(base_dir):
        sys.exit(f"ERROR: directory '{base_dir}' does not exist.")

    args.framework = args.framework or autodetect_framework(base_dir)
    args.frame = args.frame or autodetect_frame(base_dir)
    args.raspa_dir = args.raspa_dir or base_dir
    args.out = args.out or autodetect_out_path(base_dir, args.framework)

    pseudo_atoms_path = os.path.join(args.raspa_dir, "pseudo_atoms.def")
    mixing_rules_path = os.path.join(args.raspa_dir, "force_field_mixing_rules.def")
    if not os.path.isfile(pseudo_atoms_path):
        sys.exit(f"ERROR: {pseudo_atoms_path} not found.")
    if not os.path.isfile(mixing_rules_path):
        sys.exit(f"ERROR: {mixing_rules_path} not found.")

    pseudo_atoms = parse_pseudo_atoms(pseudo_atoms_path)
    mixing_rules = parse_mixing_rules(mixing_rules_path)

    mol_def_files = find_molecule_def_files(args.raspa_dir)
    print(f"Found molecule definitions: {list(mol_def_files.keys())}")

    if not args.species:
        args.species = autodetect_species_from_filename(args.frame, mol_def_files)
        if not args.species:
            sys.exit(f"ERROR: could not auto-detect the adsorbate species from frame filename "
                      f"'{args.frame}'. Pass --species explicitly (e.g. --species CO2).")

    if args.species not in mol_def_files:
        sys.exit(f"ERROR: species '{args.species}' has no matching {args.species}.def "
                  f"in --raspa-dir '{args.raspa_dir}'.")

    print(f"\n--- Using ---")
    print(f"  Framework file : {args.framework}")
    print(f"  Frame file     : {args.frame}")
    print(f"  Species        : {args.species}")
    print(f"  RASPA dir      : {args.raspa_dir}")
    print(f"  Output file    : {args.out}\n")

    fw = parse_framework_lammps(args.framework)
    hc = fw["header_counts"]
    n_atoms_fw = hc["atoms"]
    n_bonds_fw = hc.get("bonds", 0)
    n_atypes_fw = hc["atom types"]
    n_btypes_fw = hc.get("bond types", 0)

    atoms_block = fw["blocks"]["Atoms"]
    bonds_block = fw["blocks"]["Bonds"]
    assert len(atoms_block) == n_atoms_fw, f"Framework atom count mismatch: {len(atoms_block)} vs header {n_atoms_fw}"
    assert len(bonds_block) == n_bonds_fw, f"Framework bond count mismatch: {len(bonds_block)} vs header {n_bonds_fw}"

    # The frame file only ever contains one adsorbate species, so parse its
    # .def file first to learn how many atoms make up one molecule -- that's
    # what lets us regroup the frame's flat atom list into molecules.
    species_mdef = parse_molecule_def(mol_def_files[args.species])
    restart = build_restart_like_from_frame(args.frame, args.species, species_mdef["n_atoms"])
    check_box_consistency(restart["cell_vectors"], fw["box_lines"])

    # NOTE: extracted movie frames never contain framework atom positions
    # (they're per-component adsorbate movies), so framework coordinates
    # always come from the framework data file itself -- no override branch
    # needed here (unlike the old restart-based flow).

    # ------------------------------------------------------------------
    # Resolve per-component atom-type mapping and connectivity from .def files
    # ------------------------------------------------------------------
    comp_defs = {0: dict(species_mdef)}
    comp_defs[0]["name"] = args.species

    # assign new LAMMPS atom types for every unique pseudo-atom name encountered
    type_name_to_id = {}
    next_type_id = n_atypes_fw
    for comp_id, mdef in comp_defs.items():
        for name in mdef["atom_type_by_index"].values():
            if name not in type_name_to_id:
                next_type_id += 1
                type_name_to_id[name] = next_type_id

    # assign one new bond type per unique bond keyword encountered across all species
    bond_kw_to_type = {}
    next_bondtype_id = n_btypes_fw
    for comp_id, mdef in comp_defs.items():
        for (a1, a2, kw) in mdef["bonds"]:
            if kw not in bond_kw_to_type:
                next_bondtype_id += 1
                bond_kw_to_type[kw] = next_bondtype_id

    # ------------------------------------------------------------------
    # Build new Atoms / Bonds entries
    # ------------------------------------------------------------------
    new_atom_lines = []
    new_bond_lines = []
    atom_id = n_atoms_fw
    bond_id = n_bonds_fw
    mol_id = 1000  # start CO2-style adsorbate molecule IDs well above typical framework mol-IDs

    missing_params = set()

    for comp_id, mdef in comp_defs.items():
        atom_type_by_index = mdef["atom_type_by_index"]
        n_mol = restart["components"][comp_id]["n_molecules"]
        for mol in range(n_mol):
            first_id_by_local_idx = {}
            for local_idx, pname in atom_type_by_index.items():
                key = (comp_id, mol, local_idx)
                if key not in restart["adsorbate_positions"]:
                    continue
                x, y, z = restart["adsorbate_positions"][key]
                atype = type_name_to_id[pname]
                pa = pseudo_atoms.get(pname, {})
                charge = pa.get("charge", 0.0)
                atom_id += 1
                first_id_by_local_idx[local_idx] = atom_id
                new_atom_lines.append(
                    f"{atom_id:8d} {mol_id:8d} {atype:8d} {charge:10.5f} {x:12.6f} {y:12.6f} {z:12.6f}"
                )
            for (a1, a2, kw) in mdef["bonds"]:
                if a1 in first_id_by_local_idx and a2 in first_id_by_local_idx:
                    bond_id += 1
                    btype = bond_kw_to_type[kw]
                    new_bond_lines.append(f"{bond_id:8d} {btype:8d} {first_id_by_local_idx[a1]:8d} {first_id_by_local_idx[a2]:8d}")
            mol_id += 1

    n_atoms_total = n_atoms_fw + len(new_atom_lines)
    n_bonds_total = n_bonds_fw + len(new_bond_lines)
    n_atypes_total = next_type_id
    n_btypes_total = max(next_bondtype_id, n_btypes_fw)

    # ------------------------------------------------------------------
    # Masses / Pair Coeffs / Bond Coeffs additions
    # ------------------------------------------------------------------
    new_masses = list(fw["blocks"]["Masses"])
    new_paircoeffs = list(fw["blocks"]["Pair Coeffs"])
    for name, tid in sorted(type_name_to_id.items(), key=lambda kv: kv[1]):
        mass = pseudo_atoms.get(name, {}).get("mass")
        if mass is None:
            missing_params.add(name)
            mass = 0.0
        new_masses.append(f"    {tid:<4d} {mass:.6f}  # {name}")

        mr = mixing_rules.get(name)
        if mr and mr["style"] == "lennard-jones":
            eps_kcal = mr["eps_K"] * KB_KCAL
            new_paircoeffs.append(f"    {tid:<4d} {eps_kcal:.6f}  {mr['sigma']:.4f}  # {name}")
        else:
            missing_params.add(name)
            new_paircoeffs.append(f"    {tid:<4d} 0.000000  1.0000  # {name} (NO LJ PARAMS FOUND - CHECK force_field_mixing_rules.def)")

    new_bondcoeffs = list(fw["blocks"]["Bond Coeffs"])
    for kw, tid in sorted(bond_kw_to_type.items(), key=lambda kv: kv[1]):
        new_bondcoeffs.append(f"    {tid:<4d} {args.bond_k:.6f}  1.160000  # {kw} (placeholder harmonic; use fix shake / fix rigid for true rigidity, and set correct r0 for this bond)")

    if missing_params:
        print(f"WARNING: no LJ/mass parameters found for: {sorted(missing_params)}. "
              f"Check spelling against pseudo_atoms.def / force_field_mixing_rules.def.")

    # ------------------------------------------------------------------
    # Write output
    # ------------------------------------------------------------------
    with open(args.out, "w") as f:
        f.write("LAMMPS data file: framework + final GCMC-loaded adsorbate(s), generated by raspa2lammps.py\n\n")
        f.write(f"        {n_atoms_total} atoms\n")
        f.write(f"        {n_bonds_total} bonds\n")
        f.write(f"        {hc.get('angles', 0)} angles\n")
        f.write(f"        {hc.get('dihedrals', 0)} dihedrals\n")
        f.write(f"        {hc.get('impropers', 0)} impropers\n\n")
        f.write(f"           {n_atypes_total} atom types\n")
        f.write(f"           {n_btypes_total} bond types\n")
        f.write(f"           {hc.get('angle types', 0)} angle types\n")
        f.write(f"           {hc.get('dihedral types', 0)} dihedral types\n")
        f.write(f"           {hc.get('improper types', 0)} improper types\n")
        for bl in fw["box_lines"]:
            f.write(bl.strip() + "\n")

        f.write("\nMasses\n\n")
        f.write("\n".join(new_masses) + "\n")

        f.write("\nBond Coeffs\n\n")
        f.write("\n".join(new_bondcoeffs) + "\n")

        if fw["blocks"]["Angle Coeffs"]:
            f.write("\nAngle Coeffs\n\n")
            f.write("\n".join(fw["blocks"]["Angle Coeffs"]) + "\n")
        if fw["blocks"]["Dihedral Coeffs"]:
            f.write("\nDihedral Coeffs\n\n")
            f.write("\n".join(fw["blocks"]["Dihedral Coeffs"]) + "\n")
        if fw["blocks"]["Improper Coeffs"]:
            f.write("\nImproper Coeffs\n\n")
            f.write("\n".join(fw["blocks"]["Improper Coeffs"]) + "\n")

        f.write("\nPair Coeffs\n\n")
        f.write("\n".join(new_paircoeffs) + "\n")

        f.write("\nAtoms\n\n")
        f.write("\n".join(atoms_block) + "\n")
        f.write("\n".join(new_atom_lines) + "\n")

        f.write("\nBonds\n\n")
        f.write("\n".join(bonds_block) + "\n")
        f.write("\n".join(new_bond_lines) + "\n")

        if fw["blocks"]["Angles"]:
            f.write("\nAngles\n\n")
            f.write("\n".join(fw["blocks"]["Angles"]) + "\n")
        if fw["blocks"]["Dihedrals"]:
            f.write("\nDihedrals\n\n")
            f.write("\n".join(fw["blocks"]["Dihedrals"]) + "\n")
        if fw["blocks"]["Impropers"]:
            f.write("\nImpropers\n\n")
            f.write("\n".join(fw["blocks"]["Impropers"]) + "\n")

    print(f"\nWrote {args.out}")
    print(f"Total atoms: {n_atoms_total} ({n_atoms_fw} framework + {len(new_atom_lines)} adsorbate)")
    print(f"Total bonds: {n_bonds_total} ({n_bonds_fw} framework + {len(new_bond_lines)} adsorbate)")
    print(f"New atom types: {type_name_to_id}")
    print(f"New bond types: {bond_kw_to_type}")


if __name__ == "__main__":
    main()
