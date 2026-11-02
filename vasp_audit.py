"""Read-only VASP audits. No VASP license, network client or scheduler needed."""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import re
import sys
from dataclasses import dataclass
import numpy as np

NUMBER = r"[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?"
TOTEN = re.compile(r"^\s*free {2,}energy\s+TOTEN\s*=\s*(" + NUMBER + r")")
FREQUENCY = re.compile(r"^\s*(\d+)\s+f(?P<imag>/i)?\s*=.*?(" + NUMBER + r")\s+cm-1\s+(" + NUMBER + r")\s+meV")
KEYS = ("ENCUT", "PREC", "ISMEAR", "SIGMA", "IVDW", "LREAL", "LDAU",
        "LDAUTYPE", "LDAUL", "LDAUU", "LDAUJ", "ISPIN", "ISYM",
        "LDIPOL", "IDIPOL", "EDIFF")

def number(value):
    return float(value.replace("D", "E").replace("d", "e"))

def digest(path):
    if not path.is_file():
        return None
    result = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()

@dataclass
class Structure:
    cell: np.ndarray
    positions: np.ndarray
    mobile: np.ndarray
    species: list[str]

def read_poscar(path):
    lines = Path(path).read_text().splitlines()
    scales = np.array([number(x) for x in lines[1].split()])
    raw = np.array([[number(x) for x in line.split()[:3]] for line in lines[2:5]])
    if scales.size == 1:
        value = scales[0]
        if value == 0 or abs(np.linalg.det(raw)) < 1e-12:
            raise ValueError("Zero scale or singular lattice")
        scale = (-value / abs(np.linalg.det(raw))) ** (1/3) if value < 0 else value
        factors = np.repeat(scale, 3)
    elif scales.size == 3 and np.all(scales > 0):
        factors = scales
    else:
        raise ValueError("POSCAR requires one nonzero scale or three positive scales")
    cell = raw * factors
    if abs(np.linalg.det(cell)) < 1e-12:
        raise ValueError("Singular lattice")
    index = 5
    tokens = lines[index].split()
    if all(re.fullmatch(r"\d+", token) for token in tokens):
        labels = ["X" + str(i+1) for i in range(len(tokens))]
    else:
        labels = tokens
        index += 1
    counts = [int(x) for x in lines[index].split()]
    if len(counts) != len(labels) or not counts or any(n < 0 for n in counts):
        raise ValueError("Invalid species counts")
    index += 1
    selective = lines[index].strip().lower().startswith("s")
    if selective:
        index += 1
    cartesian = lines[index].strip().lower().startswith(("c", "k"))
    index += 1
    count = sum(counts)
    if len(lines) < index + count or count == 0:
        raise ValueError("Truncated or empty POSCAR")
    positions, mobile = [], []
    for line in lines[index:index+count]:
        tokens = line.split()
        positions.append([number(x) for x in tokens[:3]])
        if selective:
            flags = tokens[3:6]
            if len(flags) != 3 or any(x.upper() not in ("T", "F") for x in flags):
                raise ValueError("Invalid Selective Dynamics flags")
            mobile.append([x.upper() == "T" for x in flags])
        else:
            mobile.append([True, True, True])
    positions = np.array(positions)
    positions = positions * factors if cartesian else positions @ cell
    return Structure(cell, positions, np.array(mobile, dtype=bool),
                     [label for label, n in zip(labels, counts) for _ in range(n)])

def minimum_image(vector, cell, max_images=1000000):
    """Exact bounded closest-image search; reject excessively skewed cells."""
    vector, cell = np.asarray(vector, float), np.asarray(cell, float)
    inverse = np.linalg.inv(cell)
    fractional = vector @ inverse
    guess = fractional - np.rint(fractional)
    best = guess @ cell
    radius = np.linalg.norm(best)
    bounds = radius * np.linalg.norm(inverse, axis=0) + 1e-10
    low = np.ceil(fractional - bounds).astype(int)
    high = np.floor(fractional + bounds).astype(int)
    total = math.prod(int(b-a+1) for a, b in zip(low, high))
    if total > max_images:
        raise ValueError("Cell is too skewed for bounded image search")
    for shift in itertools.product(*(range(int(a), int(b)+1) for a, b in zip(low, high))):
        trial = (fractional - shift) @ cell
        if trial @ trial < best @ best:
            best = trial
    return best

def projected_forces(forces, mobile, cell):
    """Euclidean projection onto allowed direct-lattice directions."""
    forces = np.asarray(forces, float)
    result = np.zeros_like(forces)
    if forces.shape != mobile.shape:
        raise ValueError("Atom/force count mismatch")
    for i, flags in enumerate(mobile):
        if np.all(flags):
            result[i] = forces[i]
        elif np.any(flags):
            basis, _ = np.linalg.qr(np.asarray(cell)[flags].T, mode="reduced")
            result[i] = basis @ (basis.T @ forces[i])
    return result

def read_incar(path):
    result = {}
    if not Path(path).is_file():
        return result
    for line in Path(path).read_text().splitlines():
        line = re.split(r"[!#]", line, maxsplit=1)[0]
        for part in line.split(";"):
            if "=" in part:
                key, value = part.split("=", 1)
                result[key.strip().upper()] = " ".join(value.split()).upper()
    return result

def read_outcar(path, natoms, initial_cell):
    frames, modes, warnings = [], [], []
    electronic = "unknown"
    current_cell = np.array(initial_cell, float)
    normal_end = geometry_message = soft_stop = False
    ediff_hits = 0
    nul_count = 0
    if not Path(path).is_file():
        raise ValueError("OUTCAR is missing")
    with Path(path).open(errors="replace") as handle:
        for line in handle:
            nul_count += line.count("\0")
            if "because EDIFF is reached" in line:
                electronic = "converged"
                ediff_hits += 1
            elif "EDIFF was not reached" in line:
                electronic = "unconverged"
            normal_end |= "General timing and accounting" in line
            geometry_message |= "reached required accuracy" in line
            soft_stop |= "soft stop encountered" in line.lower()
            if "direct lattice vectors" in line and "reciprocal" in line:
                try:
                    current_cell = np.array([[number(x) for x in next(handle).split()[:3]] for _ in range(3)])
                except (StopIteration, ValueError):
                    warnings.append("Incomplete lattice block")
            if re.search(r"POSITION\s+TOTAL-FORCE", line):
                rows = []
                try:
                    while len(rows) < natoms:
                        row = next(handle)
                        if not row.strip() or set(row.strip()) == {"-"}:
                            continue
                        values = row.split()
                        if len(values) < 6:
                            break
                        rows.append([number(x) for x in values[:6]])
                except (StopIteration, ValueError):
                    pass
                if len(rows) == natoms:
                    data = np.array(rows)
                    frames.append(dict(positions=data[:, :3], forces=data[:, 3:],
                                       cell=current_cell.copy(), energy=None, scf=electronic))
                else:
                    warnings.append("Incomplete force block ignored")
            energy = TOTEN.match(line)
            if energy and frames and frames[-1]["energy"] is None:
                frames[-1]["energy"] = number(energy.group(1))
                electronic = "unknown"
            mode = FREQUENCY.match(line)
            if mode:
                if int(mode.group(1)) == 1:
                    modes = []
                modes.append(dict(index=int(mode.group(1)),
                                  imaginary=bool(mode.group("imag")),
                                  cm1=(-1 if mode.group("imag") else 1) * abs(number(mode.group(3))),
                                  energy_meV=abs(number(mode.group(4)))))
    if nul_count:
        warnings.append("NUL bytes detected; investigate concurrent writers or corruption")
    return dict(frames=frames, modes=modes, normal_end=normal_end,
                geometry_message=geometry_message, soft_stop=soft_stop,
                ediff_hits=ediff_hits, nul_count=nul_count, warnings=warnings)

def frame_metrics(frame, structure):
    projected = projected_forces(frame["forces"], structure.mobile, frame["cell"])
    norms = np.linalg.norm(projected, axis=1)
    active = np.flatnonzero(np.any(structure.mobile, axis=1))
    top = sorted(active, key=lambda i: norms[i], reverse=True)[:6]
    return dict(energy_eV=frame["energy"], scf=frame["scf"],
                fmax_projected_eVA=float(max(norms, default=0)),
                fmax_all_atoms_eVA=float(np.linalg.norm(frame["forces"], axis=1).max()),
                top_mobile=[dict(atom=int(i+1), element=structure.species[i],
                                 force_eVA=float(norms[i])) for i in top])

def load_run(directory):
    directory = Path(directory)
    structure = read_poscar(directory / "POSCAR")
    settings = read_incar(directory / "INCAR")
    output = read_outcar(directory / "OUTCAR", len(structure.positions), structure.cell)
    return directory, structure, settings, output

def audit(directory, threshold=0.03):
    directory, structure, settings, output = load_run(directory)
    ib = settings.get("IBRION", "")
    phonon = ib in ("5", "6")
    static = settings.get("NSW") == "0" or ib == "-1"
    kind = "partial_hessian" if phonon else "static" if static else "ionic"
    frames = output["frames"]
    complete = [f for f in frames if f["energy"] is not None]
    selected = complete[0] if phonon and complete else complete[-1] if complete else None
    metric = frame_metrics(selected, structure) if selected else None
    exit_file = directory / "EXIT_CODE"
    exit_code = exit_file.read_text().strip() if exit_file.is_file() else None
    modes = output["modes"]
    warnings = list(output["warnings"])
    if exit_code == "0" and any(f["scf"] != "converged" for f in complete):
        warnings.append("Exit code 0 does not establish SCF convergence")
    if phonon:
        warnings.append("Force metric is the undisplaced center, not the final displaced geometry")
    if kind == "ionic" and "2" == ib:
        warnings.append("IBRION=2 is ordinary relaxation; directory names do not establish IRC")
    return dict(kind=kind, atom_count=len(structure.positions),
                mobile_atoms=int(np.any(structure.mobile, axis=1).sum()),
                mobile_dof=int(structure.mobile.sum()), settings=settings,
                program_terminated_normally=output["normal_end"], exit_code=exit_code,
                scf_complete_force_frames=len(complete),
                scf_converged_frames=sum(f["scf"] == "converged" for f in complete),
                last_frame_scf=complete[-1]["scf"] if complete else "unknown",
                geometry_convergence_message=output["geometry_message"],
                soft_stop=output["soft_stop"],
                selected_frame=("center" if phonon else "last") if selected else None,
                force_metric=metric,
                threshold_eVA=threshold,
                selected_subspace_below_threshold=(metric["fmax_projected_eVA"] < threshold
                    if metric and metric["scf"] == "converged" and not output["nul_count"] else None),
                frequencies_cm1=[m["cm1"] for m in modes],
                imaginary_count=sum(m["imaginary"] for m in modes),
                zpe_real_modes_eV=sum(m["energy_meV"] for m in modes if not m["imaginary"])/2000 if modes else None,
                source_sha256={name: digest(directory/name) for name in ("POSCAR", "INCAR", "KPOINTS")},
                warnings=warnings)

def compare(left, right, left_frame="last", right_frame="first", tolerance=1e-5):
    ld, ls, li, lo = load_run(left)
    rd, rs, ri, ro = load_run(right)
    def select(output, which):
        frames = [f for f in output["frames"] if f["energy"] is not None]
        if not frames:
            raise ValueError("No complete force/energy frame")
        return frames[0] if which == "first" else frames[-1]
    lf, rf = select(lo, left_frame), select(ro, right_frame)
    if ls.species != rs.species:
        return dict(comparable=False, energy_difference_eV=None, issues=["Atom order/species mismatch"])
    cell_error = float(np.max(np.abs(lf["cell"] - rf["cell"])))
    displacement = max(float(np.linalg.norm(minimum_image(a-b, lf["cell"])))
                       for a, b in zip(lf["positions"], rf["positions"]))
    different = {k: [li.get(k), ri.get(k)] for k in KEYS if li.get(k) != ri.get(k)}
    unknown = [k for k in KEYS if k not in li or k not in ri]
    files = {k: [digest(ld/k), digest(rd/k)] for k in ("POTCAR", "KPOINTS")}
    issues = []
    if cell_error > tolerance or displacement > tolerance:
        issues.append("Different evaluated geometry; do not subtract endpoint energies")
    if different or unknown:
        issues.append("Numerical settings differ or are missing")
    if any(a is None or a != b for a, b in files.values()):
        issues.append("Reference file hashes differ or are unavailable")
    if lo["nul_count"] or ro["nul_count"]:
        issues.append("NUL-corrupted output cannot establish a reliable comparison")
    if lf["scf"] != "converged" or rf["scf"] != "converged":
        issues.append("SCF convergence not established")
    return dict(comparable=not issues, evaluated_frames=[left_frame, right_frame],
                max_cell_difference_A=cell_error, max_atom_difference_A=displacement,
                selected_masks_equal=bool(np.array_equal(ls.mobile, rs.mobile)),
                different_settings=different, missing_settings=unknown, reference_hashes=files,
                energies_eV=[lf["energy"], rf["energy"]],
                energy_difference_eV=rf["energy"]-lf["energy"] if not issues else None,
                spin_constraints=[li.get("NUPDOWN"), ri.get("NUPDOWN")], issues=issues)

def processes(proc="/proc", expected_ranks=None):
    groups = {}
    for entry in Path(proc).iterdir():
        if not entry.name.isdigit():
            continue
        try:
            exe = (entry/"exe").resolve().name
            if "vasp" not in exe.lower():
                continue
            cwd = os.readlink(entry/"cwd")
            match = re.search(r"^Cpus_allowed_list:\s*(.*)$",
                              (entry/"status").read_text(), re.M)
            cores = set()
            for token in match.group(1).split(",") if match else []:
                ends = [int(x) for x in token.strip().split("-")]
                cores.update(range(ends[0], ends[-1]+1))
            group = groups.setdefault(cwd, dict(cwd=cwd, pids=[], cpus=set()))
            group["pids"].append(int(entry.name))
            group["cpus"].update(cores)
        except (OSError, ValueError):
            continue
    result = []
    for cwd, group in groups.items():
        result.append(dict(cwd=cwd, rank_count=len(group["pids"]), pids=sorted(group["pids"]),
                           cpus=sorted(group["cpus"]),
                           unexpected_rank_count=(len(group["pids"]) != expected_ranks
                               if expected_ranks is not None else None)))
    overlaps = [dict(left=a["cwd"], right=b["cwd"], cpus=sorted(set(a["cpus"]) & set(b["cpus"])))
                for a, b in itertools.combinations(result, 2) if set(a["cpus"]) & set(b["cpus"])]
    return dict(jobs=result, affinity_overlaps=overlaps,
                note="Affinity overlap is a scheduling warning, not proof of CPU contention. No processes are changed.")

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("audit"); p.add_argument("directory"); p.add_argument("--threshold", type=float, default=0.03)
    p = sub.add_parser("compare"); p.add_argument("left"); p.add_argument("right")
    p.add_argument("--left-frame", choices=["first", "last"], default="last")
    p.add_argument("--right-frame", choices=["first", "last"], default="first")
    p = sub.add_parser("distance"); p.add_argument("poscar"); p.add_argument("atom1", type=int); p.add_argument("atom2", type=int)
    p = sub.add_parser("processes"); p.add_argument("--expected-ranks", type=int)
    args = parser.parse_args(argv)
    try:
        if args.command == "audit":
            result = audit(args.directory, args.threshold)
        elif args.command == "compare":
            result = compare(args.left, args.right, args.left_frame, args.right_frame)
        elif args.command == "distance":
            s = read_poscar(args.poscar)
            if min(args.atom1, args.atom2) < 1 or max(args.atom1, args.atom2) > len(s.positions):
                raise ValueError("Atom indices are 1-based and must exist")
            result = dict(distance_A=float(np.linalg.norm(minimum_image(
                s.positions[args.atom1-1]-s.positions[args.atom2-1], s.cell))))
        else:
            result = processes(expected_ranks=args.expected_ranks)
        print(json.dumps(result, indent=2, allow_nan=False))
        return 0
    except (OSError, ValueError, IndexError, np.linalg.LinAlgError) as error:
        print(json.dumps({"error": str(error)}), file=sys.stderr)
        return 2

if __name__ == "__main__":
    raise SystemExit(main())
