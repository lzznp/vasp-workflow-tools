# vasp-workflow-tools

Read-only VASP run audits and practical engineering notes from long-running surface calculations.

**中文说明**：把“程序退出、电子收敛、几何收敛”分开；按真实 Selective Dynamics 计算力；在周期晶胞中比较结构；将能量配回实际计算帧。工具不启动、不停止、不重启 VASP，也不自动读取 SSH 密钥或联网。

## Install and run

Python 3.10+ and NumPy are required.

```sh
python -m pip install .
vasp-audit audit /path/to/run --threshold 0.03
vasp-audit distance /path/to/POSCAR 1 2
vasp-audit compare /path/to/S1 /path/to/S2 --left-frame last --right-frame first
vasp-audit processes --expected-ranks 24
python -m unittest -v
```

Without installation: `python vasp_audit.py audit /path/to/run`.
All commands print JSON. `processes` requires Linux /proc; the other commands run on macOS/Linux/Windows.

## What is checked

| Command | Purpose |
|---|---|
| audit | Separate termination, electronic convergence, geometry messages and force thresholds |
| distance | Exact bounded minimum-image distance for a nonorthogonal periodic cell |
| compare | Compare evaluated OUTCAR frames, numerical inputs and reference hashes before subtracting energies |
| processes | Group real VASP processes by working directory and report overlapping CPU affinities |

- Frozen atoms cannot dominate the mobile-subspace Fmax. Mixed component flags use the **direct-lattice** directions, including Cartesian POSCARs.
- Partial Hessians report the undisplaced center's force/energy, imaginary modes and positive-real-mode ZPE.
- A zero exit code, an exhausted NSW/NELM limit, a soft stop, or a directory called TS/IRC/MECP is not scientific acceptance.
- Comparison deliberately returns no energy difference when geometry, SCF, numerical metadata or reference files are incompatible or unknown.
- OUTCAR reading is line based. Incomplete force blocks and NUL bytes produce warnings.
- This is an engineering audit, not a substitute for transition-state/endpoint validation. The projected force is a clearly defined analysis metric, not a reimplementation of VASP's stopping algorithm.

Read [工程经验与流程](ENGINEERING_NOTES.md), [输入与输出约定](FORMAT_AND_LIMITS.md), and [发布边界](PUBLICATION_SCOPE.md).

## Tests and examples

All tests build tiny synthetic cells and artificial log fragments in temporary directories. They include boundary-crossing bonds, skewed cells, frozen-atom forces, imaginary modes, unconverged exit=0 jobs, and mismatched reference frames.

No licensed VASP code, POTCAR data, private calculations, unpublished structures, experimental/patent documents, host names or credentials are distributed.

License: MIT for this repository's code and original notes. VASP remains separately licensed.
