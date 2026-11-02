# vasp-workflow-tools

A small command-line toolkit for checking VASP runs: convergence, forces on movable atoms, periodic distances, energy comparisons and CPU affinity.

中文工程笔记见 [ENGINEERING_NOTES.md](ENGINEERING_NOTES.md)，包括长时任务续算、帧与能量配对、并行调度和存储校验。

## Quick start

Requires Python 3.10+ and NumPy.

```sh
python -m pip install .

vasp-audit audit /path/to/run --threshold 0.03
vasp-audit distance /path/to/POSCAR 1 2
vasp-audit compare /path/to/S1 /path/to/S2 --left-frame last --right-frame first
vasp-audit processes --expected-ranks 24
```

You can also run the script directly: `python vasp_audit.py audit /path/to/run`.

All commands return JSON. Atom indices start at 1. Process inspection uses Linux `/proc`; the file-based commands also work on macOS and Windows. Commands only read files and process information.

## Commands

| Command | What it reports |
|---|---|
| `audit` | Program termination, SCF convergence, geometry messages, movable-subspace Fmax and frequencies |
| `distance` | Minimum-image distance in a periodic cell, including skewed cells |
| `compare` | Geometry, calculation settings and reference-file checks, followed by a TOTEN difference when comparable |
| `processes` | VASP processes grouped by working directory, rank counts and CPU affinity overlaps |

The force audit reads the actual Selective Dynamics flags. Fully frozen atoms are excluded; mixed flags use a projection onto the allowed direct-lattice directions. For partial Hessians, the audit uses the undisplaced center frame and lists imaginary modes separately from real-mode ZPE.

`compare` uses structures evaluated in OUTCAR. If geometry, SCF status or reference settings cannot be matched, it returns no energy difference and reports the reason. See [formats and limits](FORMAT_AND_LIMITS.md) for the parsing rules.

## Tests

```sh
python -m unittest -v
```

Tests use synthetic structures and log fragments. Cases cover periodic boundaries, skewed cells, frozen-atom forces, partial Hessians, incomplete output and incompatible energy references.

## License

MIT for the code and original notes. VASP and PAW datasets have their own licenses; neither is included. See [repository contents](PUBLICATION_SCOPE.md).
