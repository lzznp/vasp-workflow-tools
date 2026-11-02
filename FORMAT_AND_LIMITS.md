# Input, output and limitations

- POSCAR supports VASP 4/5 atom counts, Direct/Cartesian positions, scalar/negative-volume/three-component scales and Selective Dynamics. Velocity and predictor-corrector blocks are not changed or interpreted.
- Atom indices in the CLI are 1-based. The POSCAR species line is a comment in VASP's actual species selection; align it with your POTCAR before relying on element labels.
- Audits require an explicit POSCAR. Forces are Cartesian; mixed SD flags define a span of direct lattice directions. Euclidean projection onto that span is reported as an analysis metric.
- A complete force table is paired with the following double-space `free  energy TOTEN` summary. Single-space electronic iteration lines are excluded. No energy/force table is accepted when truncated.
- SCF status relies on explicit OUTCAR EDIFF messages. Missing messages remain unknown. Do not infer convergence from a small printed energy difference.
- Positive-real-mode ZPE is 0.5 times the sum of printed real mode energies, converted from meV to eV. Imaginary modes are identified by `f/i`, not by a minus sign in VASP's printed magnitude.
- For IBRION=5/6 the selected audit frame is the first complete, undisplaced center. Displacement convergence and the expected number of modes still require scrutiny.
- Variable-cell frames use direct lattice vectors from the OUTCAR when present. Users should verify output format and lattice association for specialized VASP modes.
- `compare` defaults to left last complete frame and right first complete frame. Geometry is periodic and atom-order preserving. No atom permutation, surface alignment, rigid rotation or chemical mapping is inferred.
- Comparison metadata intentionally requires explicit keys. Unspecified defaults return non-comparable until independently resolved; hashes are conservative and can reject scientifically equivalent differently formatted KPOINTS.
- The comparison is an electronic TOTEN difference, not a free energy, NEB barrier, spin transition probability, or MECP result.
- Linux process inspection can miss processes hidden by permission boundaries or with non-VASP executable names. An affinity overlap is a warning, not proof of resource contention. No scheduler, SSH, keychain or process control operation is performed.
- Validated with synthetic regression tests and separate local read-only smoke checks; not a comprehensive VASP parser or certified scientific validation service.
