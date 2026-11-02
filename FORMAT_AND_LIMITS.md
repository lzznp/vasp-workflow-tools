# Formats and limits

## Structures and forces

POSCAR parsing supports VASP 4/5 atom counts, Direct and Cartesian positions, scalar scaling, negative target-volume scaling, three-component scaling and Selective Dynamics. Velocity and predictor-corrector blocks are left uninterpreted.

CLI atom indices are 1-based. Element labels come from the POSCAR species line; check their order against POTCAR, which determines the species used by VASP.

An audit requires POSCAR. Forces are Cartesian, while mixed Selective Dynamics flags span direct-lattice directions, including for Cartesian input positions. The reported force metric is the Euclidean projection onto that span; it is separate from VASP's internal stopping test.

Variable-cell frames use OUTCAR lattice vectors when available. Check the lattice association when working with specialized VASP output formats.

## OUTCAR frames and frequencies

Each complete force table is paired with the following double-space `free  energy TOTEN` summary. Electronic-iteration energy lines are excluded. Incomplete tables are skipped, and NUL bytes produce warnings.

SCF convergence comes from explicit EDIFF messages. Missing messages give an unknown status, even when printed energy changes are small.

For IBRION=5/6, the audit selects the first complete, undisplaced center frame. Check displacement convergence and the expected number of modes separately. Imaginary modes are identified by `f/i`. Real-mode ZPE is half the sum of printed real-mode energies, converted from meV to eV; imaginary modes are listed separately.

## Energy comparisons

By default, `compare` takes the last complete frame on the left and the first complete frame on the right. Geometry matching is periodic and preserves atom order. It does not infer atom permutations, surface alignment or rigid rotations.

The comparison requires matching geometry, converged SCF, explicit numerical settings and reference-file hashes. Missing settings leave the result non-comparable. File hashes can also reject equivalent KPOINTS files with different formatting, so inspect the reported reason before resolving a mismatch.

The resulting difference is electronic TOTEN. Free energies, reaction barriers and spin-crossing rates require further calculations.

## Process inspection and coverage

`processes` reads Linux `/proc`. Permissions or executable naming can hide processes from the report. Overlapping CPU affinities identify a possible scheduling issue; they do not measure actual contention.

The commands do not connect over SSH or control processes. Parser coverage is checked with synthetic regression cases and local read-only smoke checks. For an unfamiliar VASP output format, inspect the selected frames and warnings before using the results.
