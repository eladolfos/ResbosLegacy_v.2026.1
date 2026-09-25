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
`[grids] experimental` for "ExpCustomGrid" mode (below), `[legacy]` overrides (`bmax`, `nonpert`, `ibeam`,
`fract_n` — the latter two for fixed-target/nuclear-target runs, e.g. `ibeam = 0` + `fract_n = 0.54d0` for a
p+Cu target; applied to both the `legacy_Y` and `legacy_main` `.in` files), `[resbos]` (VEGAS settings, seed,
luminosity, output format and named cut sets `runs = atlas, nocuts` with `[cuts.<name>]`) and `[environment]`
(modules, LHAPDF, HOPPET, ROOT paths used inside the generated job scripts — set these if your installs
differ from the defaults).

### `[campaign] compute`: resNLO, NLO, or both

The driver can build two different physics results from the same grid/campaign; pick one or
both with `[campaign] compute` (default `resNLO`, so every existing `.ini` without this key
behaves exactly as before):

```ini
[campaign]
compute = resNLO           # default -- the CSS-resummed result (section 3.3's job chain)
# compute = NLO             # fixed order, via Legacy's own LTO=1/2/3 pieces + resbos + hadd
# compute = NLO, resNLO     # both in one campaign; legacy_Y (LTO=3) is computed once and shared
```

See `E201_e605_NLO_resNLO.ini` and `7TeV_WpWm_NLO_resNLO.ini` for worked examples of
`compute = NLO, resNLO` (copies of `E201_e605_fine.ini`/`7TeV_WpWm_example.ini` with that one
key added, plus the reuse note below).

#### resNLO (the default workflow)

This is the chain already described in section 3.3: Legacy's `LTO=3` (Y-piece) combines with
`w_pert`/`w_asym`'s fixed-order pieces in `get_yk_new` into a K-factor ("Yk"/`ykR`) grid, which
`resbos` reads alongside Legacy's `LTO=0` ("Main data grid" — the CSS-resummed cross section
itself) to produce weighted events. The physics content is the full CSS-resummed qT spectrum,
matched to fixed order at large qT via the K-factor — this is "resNLO" (resummed-and-matched
NLO) in the sense `w321`/Yao's naming uses it, and the result this repository has validated end
to end against real data (E605, see the E605 campaign `.ini` files' own comments/README history).
Pieces needed: `w_pert`, `w_asym`, Legacy `LTO=3` and `LTO=0`, `get_yk_new`, one `resbos` run per
`[resbos] runs` entry.

#### NLO (fixed order, via phase-space slicing)

A genuine **fixed-order** (no resummation) NLO cross section can be built from Legacy alone,
without `w_pert`/`w_asym`/`get_yk_new` at all, using the standard two-cutoff phase-space-slicing
method:

- Legacy's Y-piece is *defined* as `Y = FixedOrder(2→3, real emission) − Asymptotic` (the CSS
  formula's own large-qT expansion) — so by construction `Asymptotic (LTO=2) + Y-piece (LTO=3)
  = FixedOrder(2→3)`, valid down to some small `qT_Sep` where the asymptotic expansion is an
  accurate stand-in for the real matrix element.
- Below `qT_Sep`, the 2→3 real-emission cross section develops a soft/collinear singularity as
  qT→0; Legacy's `LTO=1` ("DeltaSigma"/`NLO_Sig`) is the regulated singular 2→2-kinematics piece
  that captures exactly this region (no Y-piece needed there — there is no "large-qT expansion"
  below the cutoff, the calculation is intrinsically 2→2).
- Run `resbos` **twice**: once with Main grid = `LTO=2` (Asymptotic) and Y grid = `LTO=3`
  (Y-piece) — the `qT>qT_Sep` piece — and once with Main grid = `LTO=1` (DeltaSigma) and **no**
  Y grid at all (the `.in` file's "Y piece grid" field set to `-`) — the `qT<qT_Sep` piece. Add
  the two resulting event samples with ROOT's `hadd`. The artificial `qT_Sep` dependence
  introduced by the slicing is designed to cancel between the two pieces in their sum (standard
  two-cutoff method), so the combined `.root` is a genuine fixed-order NLO result, independent
  of `qT_Sep` up to residual power corrections.

This is **not** boson-specific (W/Z/A0 all work the same way): `resbos_root.f` figures out which
physics mode a "Main data grid" file holds by reading an `LTOpt` field from the grid file's *own*
header (`YUAN_MAIN`, `Read(2,*) ECMC, iBeam, nDummy, LTOpt, iProc`) and sets its internal event
kinematics from that, not from `Type_V`/`iProc` alone — the boson identity is already baked into
the grid's own numbers regardless of which LTO produced it. This is exactly the mechanism that
lets `resbos` accept a `LTO=1`/`LTO=2` grid as a "Main data grid" at all (normally that field is
always `LTO=0`), and it is what the "run resbos twice + hadd" recipe above relies on.

Pieces needed for Fixed Order (no `w_pert`/`w_asym`/`get_yk_new`):

| Piece | Legacy `LTO` | Role |
|---|---|---|
| DeltaSigma (`legacy_dsi`) | `1` | Singular 2→2 piece, `qT < qT_Sep` |
| Asymptotic (`legacy_asy`) | `2` | Large-qT expansion, paired with the Y-piece for `qT > qT_Sep` |
| Y-piece (`legacy_Y`) | `3` | `FixedOrder(2→3) − Asymptotic`; **shared with resNLO** if both are requested |
| `resbos` (asy+Y run) | — | Main=`legacy_asy`, Y=`legacy_Y` → `qT>qT_Sep` events |
| `resbos` (dsi-alone run) | — | Main=`legacy_dsi`, Y=`-` → `qT<qT_Sep` events |
| `hadd` | — | Combines the two `.root` files into `<name>_<tag>_<run>_NLO.root` |

`legacy_dsi`/`legacy_asy` reuse the `legacy_Y` template file by default (`patch_stage()`
overwrites the `LTO` token regardless of what the template starts with) — set `[templates]
legacy_dsi`/`legacy_asy` only if you need a genuinely different template. `[grids] experimental`
(ExpCustomGrid mode) does not support `compute = NLO`: fixed order needs a real, gap-free
`(Q,qT,y)` rectangle the same way resNLO's `resbos` step does (see the "ExpCustomGrid" WARNING below).

The asy+Y `resbos` run also needs `iYGrid=1` in its `.in` file: `resbos_root.f` picks the Y-grid
*format* it expects from an `iYGrid` field in `resbos.in` (not by inspecting the file), and every
`resbos.in` template ships with `iYGrid=2` — the format `get_yk_new` writes (17 header lines, an
extra `R_Ai` column), tuned for resNLO's own run. Legacy's raw `LTO=3` output is the *older*
10-column/15-header layout instead (`iYGrid=1`); feeding it to a `resbos.in` still set to
`iYGrid=2` makes `resbos_root.f` print `Must have 17 comment lines in the Y-Grid file` and exit
immediately (**exit code 0, no error, no events** — easy to miss). `write_resbos_run()` sets
`iYGrid=1` automatically whenever the Y grid comes from `legacy/` instead of `get_yk_new/` (i.e.
only for NLO's asy+Y run), so this is handled for you — mentioned here in case you ever build a
`resbos.in` by hand for this workflow.

**Caveat**: the `header_lines`/`lines_per_point` shape used for `LTO=1`/`LTO=2` (15 header lines,
1 data line/point, same as `LTO=0`) was derived by reading `legacy_final_vesion/main.for`'s
write statements, not by running Legacy (no `gfortran` on this machine) — spot-check a real
`LTO=1`/`LTO=2` output on the HPCC (line counts, `merge_*_array.sb`'s row-count check) before
trusting a large production run.

### Reusing outputs already in `dest/` (don't repeat what's already computed)

Every stage (`w_pert`, `w_asym`, each Legacy `LTO`, `get_yk_new`, each `resbos` run) checks, before
generating its shard/submit lines, whether its own output already exists in `dest/` and is valid
(a Legacy/`w_pert`/`w_asym`/`get_yk_new` `.out` passes the same row-count check its merge job would
run; a `resbos` run counts as done if `dest/<job>.root` exists). If so, it's printed as "already
exists -- reusing, not resubmitted" and left out of `submit_all.sh` entirely — only stages that are
missing or whose upstream inputs actually changed get (re)submitted, and downstream
`--dependency=afterok:` clauses only reference the jobs actually queued this run. This is what
makes it cheap to layer `compute = NLO` onto a `dest/` that already has a finished `resNLO` run (or
vice versa): rerun with `--reuse` and only the genuinely new pieces are queued, `legacy_Y` included
(shared, so it's computed once total, not once per `compute` value) — see the examples above.

### Building a grid for a fixed-target campaign

`make_grid_from_data.py` has two modes for building `q_grid.inp`/`y_grid.inp` (`qt_grid.inp` is left
untouched in both — point `[templates]` at a source whose `qt_grid.inp` is already the one you want):

- `unique` — the grid points are the table's own unique Q/y values (see ExpCustomGrid below for when
  this fits vs. wastes compute).
- `generate` — a genuine fine, dense rectangular grid (N linearly- or log-spaced points over a
  continuous range) — the same role as `templates/7TeV_WpWm`'s 80 Q × 143 y grid, needed to get a real
  `resbos`/resNLO result (`resbos` interpolates over the grid; it needs it dense and complete, not just
  wherever an experiment happened to measure). The range is either typed by hand (`--q-min`/`--q-max`)
  or taken from a data table's own min/max (`--q-from TABLE`, borrowing only the table's *span*, not its
  individual values, unlike `unique` mode):
  ```bash
  python3 make_grid_from_data.py generate \
      templates/E605_fine/legacy/inp/q_grid.inp templates/E605_fine/legacy/inp/y_grid.inp \
      --n-q 40 --n-y 25 --q-min 4.0 --q-max 20.0 --y-min -0.5 --y-max 0.5
  ```
  As with `unique` mode, copy the same two output files into `w_pert`/`w_asym`'s `inp/` too so all three
  match by md5 (`run_process.py`'s grid-lockstep check).

### ExpCustomGrid: one job per experimental point

A fixed-target campaign whose measured points don't form a dense rectangle (most of them: each `(y,Q)` row
is essentially unique) wastes enormous compute under the normal one-shared-grid model — `make_grid_from_data.py`
prints how bad it would be (e.g. E866f: 15 points but 15 unique Q × 15 unique y = 225 grid points, 93% wasted;
E605 is the exception, 119 points → 18×7=126, only 6% wasted, so it uses a normal grid built with that script
instead). For the bad cases, set `[grids] experimental = /path/to/data_table` (a whitespace-separated table,
one measured point per row, `exp_y_col`/`exp_q_col` pick the columns — default 1/2, matching Yao's tables
under `New_kFactorCT25/FixedTarget_pp830a016_yao_09182026/Workspace/*/<experiment>/<experiment>`). This
switches `run_process.py` into a mode that mirrors those experiments' own `run.sh`/`combine.sh`: one job per
measured `(y,Q)` point (single-value `q_grid.inp`/`y_grid.inp`, the template's full `qt_grid.inp` swept in
every job), run as a SLURM array (`dest/<folder>/run_<job>_points.sb`, optionally throttled with `[grids]
throttle`), then merged (`run_<job>_merge.sb`, header-stripped concatenation + the same row-count check as a
normal sharded merge) into `<job>_combined.out`. `get_yk_new` and `resbos` then run once per boson, fed those
merged files (`--resbos-only` still refuses to run against an experimental-grid `dest`, since it assumes
plain, non-`_combined` file names). `get_yk_new`'s dummy `R_Ai` (NLO only) is written as the *full*
unique-Q × unique-y rectangular closure, not the real (irregular) point set — `get_yk_new.f`'s `CheckSum`
requires `NQ×Ny×NqT == total rows`, which a genuine sparse point set can never satisfy even with an
all-1 R_Ai the values of which are unused at NLO; this closure only has to pass that structural check,
it adds no physics and the K-factor itself is still computed from the real points only.

**`resbos` itself, though, generally does NOT reach a real result in this mode.** It runs its own
`CheckSum` (`resbos_root.f`, same block-counting logic as `get_yk_new.f`'s) on the *Main data grid* —
Legacy's real, non-dummy `LTO=0` output — with the same `NQ×Ny×NqT==rows` requirement. Unlike `R_Ai`,
those are real physics values `resbos` actually samples from, so they can't be padded with fake rows the
way `R_Ai` can; `run_process.py` prints a `WARNING` at prepare time when the points aren't a perfect
`Q × y` rectangle, since this only ever surfaces (`File checksum BAD. The grid file is corrupt!`) after
the full run finishes on the cluster. For a real `resbos`/`.root` result, build a genuine rectangular
grid instead (`make_grid_from_data.py` + a normal `[templates]` campaign — e.g. `E201_e605_full.ini`).

The full grid is 80 Q x 153 qT x 143 y = 1,750,320 points.

### Timing log and the final .root file

`dest/timing.log` tracks the campaign end to end: `submit_all.sh` writes the first line (`T0 <epoch>`,
the campaign's own start time) before submitting anything, and every merge job, `get_yk_new` job and
`resbos` job appends one locked (`flock`) line when it finishes (so the file only ever grows, and a
stage showing up in it means that stage is done) with two numbers:
- `elapsed_since_start` — wall-clock time since the whole campaign began (`T0`), queue wait included.
- `duration` — that stage's own time in isolation, excluding time spent waiting on earlier stages:
  for `w_pert`/`w_asym`/`legacy_Y`/`legacy_main` (array + merge), `submit_all.sh` passes down the
  moment the array was submitted (the array itself has no dependency, so that's effectively its start);
  for `get_yk_new`/`resbos` (single jobs gated by `--dependency=afterok:...` on earlier stages), it's
  just bash's own `$SECONDS`, which only starts counting once SLURM actually dispatches that job --
  i.e. after its dependency is satisfied, so the wait is already excluded for free.

Concurrent array-task output isn't logged (119 lines per stage would be too noisy — only the merge
step per stage, and `get_yk_new`/`resbos`, log a line), so the file stays a short, readable summary
you can `cat` for a quick "how far did it get, how long did each step take" check.

`resbos_root` writes its ntuple as `<jobname>.root` inside `dest/resbos/`; the `resbos` job script moves
it up to `dest/<jobname>.root` on success, so the final result sits directly in the campaign's top-level
folder — easy to find and `scp`/copy out without digging into `resbos/`.

### 3.2 Prepare the run folder

```bash
python3 run_process.py 7TeV_WpWm_example.ini
```

This creates `run_7TeV_WpWm/` with everything a run needs (the executables from `bin/`, grids, `.in` files rewritten for the
chosen energy / PDF / boson, shard scripts, SLURM scripts) and `submit_all.sh`. It prints which executable it took
from where, the grid size, and the list of jobs per boson. Nothing is submitted yet, so you can inspect the `.in`
files first. The `.in` files are generated from the W+ templates: only ECM, `lha_<pdf>`, `Type_V`/`JWTYPE`, the
active range and the options set in the `.ini` are changed.

Add `--clean-shards` to have each `w_pert`/`w_asym`/legacy merge job delete its own `<job>_shards/` folder (the
per-shard `.in`/`.out`/logs; several GB for the full grid) right after the merge's row-count check passes, freeing
the space as soon as the merged `.out` is verified good. If the check fails the shards are kept, so you can inspect
the bad shard. Without the flag the shards are kept in every case.

```bash
python3 run_process.py 7TeV_WpWm_example.ini --clean-shards
```

### 3.3 Submit

```bash
python3 run_process.py 7TeV_WpWm_example.ini --reuse --submit     # on the HPCC login node
# or, if the folder is already prepared:
bash run_7TeV_WpWm/submit_all.sh
```

(`--submit` on a fresh `dest` works without `--reuse`; use `--reuse` when the folder already exists. `--clean-shards`,
`--reuse` and `--submit` are independent flags and combine freely, e.g. `--clean-shards --submit` to prepare, submit
and have the shards cleaned up as each merge job finishes.)

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

### 3.5 Add a resbos run (e.g. a different set of cuts) without recomputing everything

`w_pert`, `w_asym`, Legacy and `get_yk_new` do not depend on `[resbos]` / `[cuts.<name>]` at all, so the best way to
get several cut sets is to list them all in `[resbos] runs` **before** the first `python3 run_process.py ... --submit`
(section 3.1) — they share the same upstream jobs, run in the same submission, and cost nothing extra upstream:

```ini
[resbos]
runs = atlas, nocuts        # one resbos job per name, both depending on the same Yk grid / Legacy main grid
```

If you only realize you need another cut set **after** a campaign has already finished (`legacy/*_main.out` and
`get_yk_new/*.out` already exist in `dest/`), don't rerun `run_process.py` normally — `submit_all.sh` always
resubmits the whole chain, wasting the ~16 h of Legacy and ~3 h of `w_pert`. Instead add the new run to the `.ini`
and use `--resbos-only`, which only (re)writes and submits the `resbos` job(s), reusing the existing upstream outputs
directly (no dependency wait, since they are already there):

```ini
[resbos]
runs = nocuts                     # add the new run (drop or keep 'atlas'; already-existing resbos/*.in are just rewritten)

[cuts.nocuts]
lepton = 0.0, -10, 0.0, 10000.0, 99.0
mt_met = 0.0, 10000., 0.0
```

```bash
python3 run_process.py 7TeV_WpWm_example.ini --resbos-only --submit
```

It re-copies the current `resbos_root` from `bin/` (so a rebuilt `resbos` is picked up too), checks that
`legacy/<...>_main.out` and `get_yk_new/<...>.out` already exist for each process in `[campaign] processes` (and
stops with a clear error if not — run the full campaign first), then writes and submits only the new `resbos/*.in` +
`run_*.sb`. Without `--submit` it just prepares them and prints the `sbatch` command to run manually.

### 3.6 Change the setup and run again

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

`W+`, `W-` and `A0` (fixed-target photon Drell-Yan, e.g. E605) are wired up end to end in the driver;
`Z0` is defined in `PROCS` but `get_yk_new.f` still `STOP`s for it (no `convert` factor implemented
there yet).

## 4. Repository layout

```
run_process.py            driver (prepares a run and writes submit_all.sh)
7TeV_WpWm_example.ini     example configuration (compute = resNLO, the default)
7TeV_WpWm_NLO_resNLO.ini  same campaign, compute = NLO, resNLO (fixed order + resummed, see section 3.1)
E201_e605_fine.ini        E605 fixed-target example (compute = resNLO)
E201_e605_NLO_resNLO.ini  same campaign, compute = NLO, resNLO
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
