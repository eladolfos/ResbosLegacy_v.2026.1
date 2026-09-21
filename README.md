# ResbosLegacy_v.2026.1

All-in-one ResBos / Legacy workflow: resummed transverse-momentum (q<sub>T</sub>) distributions of
W/Z/Higgs bosons in hadron collisions (CSS formalism, Legacy/ResBos family from MSU).
This repository holds the five codes of the pipeline, a one-file driver that prepares and chains a whole
production run on SLURM, and a ready-to-use example (W<sup>+</sup> and W<sup>-</sup> at 7 TeV).

| Folder | Executable | Role |
|---|---|---|
| `w_pert_08112022/` | `w_pert` | Fixed-order perturbative (Arnold–Reno) `dσ/dy/dqT²` on a (Q, qT, y) grid |
| `w_asym_08112022/` | `w_asym` | Asymptotic / small-qT expansion of the same, on the same grid |
| `legacy_final_vesion/` | `main` | Legacy: CSS resummed piece + Y piece + non-perturbative Sudakov, steered by a `.in` file |
| `get_yk_new/` | `get_yk_new` | Combines `w_pert`, `w_asym` and Legacy (Y piece) into the Y / `ykR` K-factor grid |
| `resbos/` | `resbos_root` | Monte Carlo event generator; reads the Legacy and Yk grids |

Data flow, per boson (W+ and W- share nothing):

```
w_pert  ------------------\
w_asym  --------------------> get_yk_new --> Yk grid --\
legacy (LTO=3, Y piece) ---/                            +--> resbos  (MC events / ntuple)
legacy (LTO=0, main grid) -----------------------------/
```

Legacy therefore runs twice per boson (the two `.in` files differ only in `LTO`).
`w_pert`, `w_asym` and Legacy must run on the **same** (Q, qT, y) grids; the driver checks this.

---

## 1. Requirements

Everything is meant to be built and run on the **MSU HPCC** (SLURM). The codes cannot be built on a machine
without these:

- `gfortran`, `g++`, `make` (the HPCC modules `GCC/13.2.0 OpenMPI/4.1.6-GCC-13.2.0 powertools` provide them)
- [LHAPDF](https://lhapdf.hepforge.org/) with `lhapdf-config` and the PDF set you want (the example uses `CT25NNLO`)
- [HOPPET](https://hoppet.hepforge.org/) (library `libhoppet_v1`), needed by Legacy
- [ROOT](https://root.cern/) with `root-config`, needed by `resbos`
- `python3` (3.6+) — **standard library only**, no conda/pip environment is needed
- SLURM (`sbatch`, `squeue`)

The default install paths used by the scripts are the author's HPCC ones:

```
LHAPDF     /mnt/home/lopezels/InstallSources/LHAPDF
HOPPET     /mnt/home/lopezels/InstallSources/HOPPET1
ROOT_SETUP /mnt/home/lopezels/InstallSources/ROOT/bin/thisroot.sh
```

If yours differ, override them as shown below (both in the build step and in the `.ini`).

## 2. Setup

### 2.1 Get the code

```bash
git clone https://github.com/eladolfos/ResbosLegacy_v.2026.1.git
cd ResbosLegacy_v.2026.1
```

### 2.2 Build the executables (`setup_resbos_legacy.sb`)

The executables that come inside `templates/` are old snapshots and are **not** used. Build fresh ones from the current sources.
`setup_resbos_legacy.sb` is plain bash (no python) and must be run **from the top folder of the repo**:

```bash
sbatch setup_resbos_legacy.sb            # compile on a compute node (2 h limit, 4 cores)
# or, directly on a login node:
bash -l setup_resbos_legacy.sb
```

What it does, for each of `w_pert`, `w_asym`, `legacy`, `get_yk_new`, `resbos`:

1. copies the top-level source files of the code's folder into `build/<code>/` (the source folders are never touched);
2. runs `module purge; module load ...`, sets the LHAPDF/HOPPET paths (and sources ROOT for `resbos`);
3. runs a from-scratch `make` there (Legacy is linked against `$HOPPET`, overriding the path hardcoded in its Makefile);
4. copies the result to `bin/<code>/<exe>`.

Progress goes to `slurm_setup_<jobid>.out`; the compiler output of each code is in `build/<code>/build.log`.
When it finishes you should have:

```
bin/w_pert/w_pert   bin/w_asym/w_asym   bin/legacy/main   bin/get_yk_new/get_yk_new   bin/resbos/resbos_root
```

Options are environment variables placed before the command:

```bash
ONLY="legacy resbos" sbatch setup_resbos_legacy.sb             # build only some codes
OUT=/some/folder sbatch setup_resbos_legacy.sb                 # write them to /some/folder/<code>/<exe> instead of bin/
                                                               # (then set [executables] dir = /some/folder in the .ini)
LHAPDF=/my/LHAPDF HOPPET=/my/HOPPET ROOT_SETUP=/my/ROOT/bin/thisroot.sh \
MODULES="GCC/13.2.0 powertools" sbatch setup_resbos_legacy.sb  # different installs / modules
NO_ENV=1 bash setup_resbos_legacy.sb                           # do not touch modules/paths: use the current shell
```

If a code fails, the script prints the last lines of its log, keeps building the others and exits with an error.

### 2.3 How the driver finds them

`run_process.py` uses `bin/<code>/<exe>` (the `bin/` next to `run_process.py`) **by default**, so after section 2.2
there is nothing to configure. The old executables stored inside `templates/` are never used. If any executable is
missing, the driver lists which ones and stops before copying anything; build them first.

To use executables from somewhere else, add to the `.ini`:

```ini
[executables]
dir = /path/to/folder              # <dir>/<code>/<exe>, the same layout as bin/ (relative paths are relative to the .ini)
# w_pert = /path/to/w_pert         # or one key per code (w_pert, w_asym, legacy, get_yk_new, resbos); a key wins over dir
```

Every run prints the exact path each executable was copied from.

### 2.4 Check the PDF set

The `.ini` names an LHAPDF set (`pdf = CT25NNLO`). It must be installed on the HPCC:

```bash
ls $(lhapdf-config --datadir) | grep CT25NNLO
```

## 3. Running the example: W+ and W- at 7 TeV

Files involved:

| File | Purpose |
|---|---|
| `7TeV_WpWm_example.ini` | The whole run in one file: energy, PDF, bosons, order, shards, optional resbos cuts |
| `run_process.py` | The driver: prepares a run folder and writes the job chain |
| `templates/7TeV_WpWm/` | The *source*: `get_yk_new/ legacy/ resbos/ w_asym/ w_pert/` with the W+ `.in` files and the `inp/` grids (its old executables are not used) |
| `shrds_scripts/` | Scripts that split the long runs along Q into SLURM job arrays and merge them back |

### 3.1 Edit the `.ini`

Open `7TeV_WpWm_example.ini`. The main settings (relative paths are relative to the `.ini`):

```ini
[paths]
source = templates/7TeV_WpWm    # template folder
dest   = run_7TeV_WpWm          # NEW folder, created by the driver; the jobs run here

[campaign]
name      = 7TeV                # prefix of every file (keep it short: the Fortran codes truncate paths at 100 chars)
ecm       = 7000                # GeV
processes = W+, W-              # each one runs the whole chain
pdf       = CT25NNLO            # LHAPDF set
order     = NLO                 # get_yk_new order; NNLO needs real R_Ai tables ([get_yk_new] r_ai_Wp / r_ai_Wm)

[executables]                   # optional; default: bin/<code>/<exe> from setup_resbos_legacy.sb (section 2.2)

[shards]                        # split along Q; 1 = one array task; resbos is never sharded
legacy = 20                     # ~16 h serial for the full grid
w_pert = 4                      # ~3 h serial
w_asym = 4
```

Optional sections (all documented in the example file): `[grids] active` to run a sub-range of the grid,
`[legacy]` overrides (`bmax`, `nonpert`), `[resbos]` (VEGAS settings, seed, luminosity, output format and named
cut sets `runs = atlas, nocuts` with `[cuts.<name>]`) and `[environment]` (modules, LHAPDF, HOPPET, ROOT paths
used inside the generated job scripts — set these if your installs differ from the defaults).

The full grid is 80 Q x 153 qT x 143 y = 1,750,320 points.

### 3.2 Prepare the run folder

```bash
python3 run_process.py 7TeV_WpWm_example.ini
```

This creates `run_7TeV_WpWm/` with everything a run needs (the executables from `bin/`, grids, `.in` files rewritten for the
chosen energy / PDF / boson, shard scripts, SLURM scripts) and `submit_all.sh`. It prints which executable it took
from where, the grid size, and the list of jobs per boson. Nothing is submitted yet, so you can inspect the `.in`
files first. The `.in` files are generated from the W+ templates: only ECM, `lha_<pdf>`, `Type_V`/`JWTYPE`, the
active range and the options set in the `.ini` are changed.

### 3.3 Submit

```bash
python3 run_process.py 7TeV_WpWm_example.ini --reuse --submit     # on the HPCC login node
# or, if the folder is already prepared:
bash run_7TeV_WpWm/submit_all.sh
```

(`--submit` on a fresh `dest` works without `--reuse`; use `--reuse` when the folder already exists.)

`submit_all.sh` submits, for each boson, the chain below with `--dependency=afterok`, so every step starts only when
the previous ones finished successfully:

```
w_pert (4 shards)  ─ merge ─┐
w_asym (4 shards)  ─ merge ─┼─> get_yk_new ─┐
legacy Y, LTO=3 (20 shards) ─ merge ─┘      ├─> resbos
legacy main, LTO=0 (20 shards) ─ merge ─────┘
```

Each merge job also checks the number of lines of the merged file and fails if it is wrong
(3 header + 1 line/point for `w_pert`/`w_asym`; 15 header + 3 lines/point for Legacy Y; 15 header + 1 line/point for
Legacy main). Jobs that need earlier outputs link them (relative symlinks; `[paths] stage_outputs = copy` makes real
copies) when they start and stop with a clear error if one is missing.

### 3.4 Monitor and collect results

```bash
squeue -u $USER
scancel <jobid>                     # to stop one
```

Logs are `slurm_*.out/.err` inside each code's folder (shard logs are in `<job>_shards/slurm_logs/`). Results:

| Step | Output (inside `run_7TeV_WpWm/`) |
|---|---|
| `w_pert`, `w_asym` | `w_pert/7TeV_w_pert_Wp.out`, `w_asym/7TeV_w_asym_Wp.out` (and `..._Wm`) |
| Legacy | `legacy/7TeV_legacy_Wp_Y.out` (LTO=3) and `legacy/7TeV_legacy_Wp_main.out` (LTO=0) |
| `get_yk_new` | `get_yk_new/7TeV_Yk_Wp_NLO.out` |
| `resbos` | in `resbos/`; see `resbos_output_<job>_<jobid>.log` for the run summary and the output files (a ROOT ntuple for `output = ROOTNT2`) |

For reference, a full-size run takes about 3 h (`w_pert`) and about 16 h serial for Legacy (divided by the number of
shards); `resbos` takes minutes for tens of millions of events.

### 3.5 Change the setup and run again

Copy the `.ini`, edit it and run the driver with the new file. For example to change the energy, PDF, order or the
resbos cuts edit `ecm`, `pdf`, `order` or the `[resbos]` / `[cuts.<name>]` sections, and use a different `dest`
(or `name`) so that you do not overwrite a previous run.

```bash
cp 7TeV_WpWm_example.ini 8TeV_WpWm.ini      # edit name, ecm, dest, ...
python3 run_process.py 8TeV_WpWm.ini --submit
```

`--reuse` refreshes an existing `dest`, but it refuses to reuse a `<job>_shards/` folder whose `.in` changed
(the shard scripts would silently keep the old energy); use a new `dest` in that case. The driver also refuses
job names that would be longer than 95 characters and grid sets that are not identical in `w_pert`, `w_asym`
and `legacy`.

Only `W+` and `W-` are wired up in the driver.

## 4. Repository layout

```
run_process.py            driver (prepares a run and writes submit_all.sh)
7TeV_WpWm_example.ini     example configuration
setup_resbos_legacy.sb    builds the five executables into bin/
shrds_scripts/            Q-sharding scripts (make_shards / merge_shards / run_*_array / submit_*_array)
templates/7TeV_WpWm/      source template for the example (.in files, grids; its old executables are not used)
w_pert_08112022/  w_asym_08112022/  legacy_final_vesion/  get_yk_new/  resbos/
                          source code of the five programs (Fortran 77 fixed format + C/C++ bridges)
bin/  build/  run_*/      generated (git-ignored): bin/ = fresh executables, run_*/ = prepared runs
```

Notes on the source: the Fortran needs `-fno-automatic` (already in every Makefile; never drop it), and
`README.md` files inside `legacy_final_vesion/` and `w_asym_08112022/` describe the individual codes and their input
files in detail. The five code folders were previously independent git repositories (still available, unchanged, under
`github.com/eladolfos/`); they are now part of this single repository.
