# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this directory is

`LegacyCode/` is **not itself a git repository**. It is a workspace holding five independent sibling git repos (each with its own `origin` under `github.com/eladolfos/`) that form one physics pipeline: resummed transverse-momentum (qT) distributions of W/Z/Higgs bosons in hadron collisions (CSS formalism, Legacy/ResBos family from MSU). Run `git` commands inside the specific subproject, never at this level.

| Dir | Executable | Role |
|---|---|---|
| `legacy_final_vesion/` (sic) | `main` | Legacy: CSS resummed piece + Y-piece + non-perturbative Sudakov; steered by `legacy.in`. `LTO=3` writes the Y-piece grid consumed downstream as `legacy_y_*.out` |
| `w_pert_08112022/` | `w_pert` | Fixed-order perturbative (Arnold–Reno) `dσ/dy/dqT²` on a (Q,qT,y) grid |
| `w_asym_08112022/` | `w_asym` | Asymptotic / small-qT expansion of the same, on the same grid |
| `get_yk_new/` | `get_yk_new` | Combines the three outputs above (+ optional `R_Ai` angular-correction table) into the Y/`ykR` K-factor grid |
| `resbos/` | `resbos_root` | Monte Carlo event generator; reads the grids produced by Legacy/`get_yk_new` via `resbos.in` |

Data flow: `w_pert` + `w_asym` + Legacy (`LTO=3`) → `get_yk_new` → `ykR` grid → `resbos`. `resbos` also reads a second Legacy run, `LTO=0` (its "Main data grid"), so Legacy runs **twice per process** (the two `.in` files differ only in `LTO`). Per process (W+, W−) the inputs differ only in `Type_V`, `JWTYPE` (+1/−1 in `w_pert`/`w_asym`) and the `resbos` grid paths; `resbos` infers the boson from the grids' `TYPE_V` header, so W+ and W− cannot share any result.

**`Example_E211_ATLAS8TeV_WpWm/`** is the reference for the current, correct workflow (ATLAS 8 TeV W±, CT25NNLO, 80×153×143 grid). Its outputs have fixed line counts, used for validation: `w_pert`/`w_asym` 3 header + 1 line/point; Legacy `LTO=0` 15 + 1 line/point; Legacy `LTO=3` 15 + 3 lines/point; `get_yk_new` output 17 + 3 lines/point. Reference timings: `w_pert` 3 h 10 min; `resbos` ~4 min for 25 M events (so Legacy, ~16 h serial, dominates). `get_yk_new/E211/` is an older, partially failed snapshot of the same campaign (its `E211_Yk_Wp_NLO.out` is a 17-line stub) — prefer the `Example_` folder.

**Subprojects that already have their own `CLAUDE.md`** (read these first when working in them): `get_yk_new/`, `resbos/`, `w_pert_08112022/`. `legacy_final_vesion/` and `w_asym_08112022/` have only a `README.md` (detailed install/run/HPCC instructions — the input-file format tables in them are accurate and worth reading before editing `.in` files).

## Environment

- Fortran 77 (fixed format, `.f`/`.for`) plus C++ (`lhapdf.cpp` bridge, `angular.cc`, `froot.c`). Everything relies on `-fno-automatic` (static locals) — never drop it.
- The current dev machine has **no `gfortran`, `lhapdf-config`, or ROOT installed** (only `g++`), so the codes cannot be built or run locally. Builds and runs happen on MSU HPCC (SLURM; paths like `/mnt/home/lopezels/InstallSources/{LHAPDF,HOPPET1}`), and outputs are copied back here.
- `legacy_final_vesion/Makefile` hardcodes `HOPPETLIBS = -L/mnt/home/lopezels/InstallSources/HOPPET1/lib -lhoppet_v1`; `hoppet/` is a bundled copy (a plain directory here, not a `.gitmodules` submodule) for building from source instead.
- Commands (per subproject, from its own directory): `make` / `make clean`. For `w_pert`, `w_asym`, Legacy the Makefile needs `lhapdf-config` on `PATH`.

## Running

Legacy, `w_pert`, and `w_asym` take a jobname argument or stdin; the input file's grid paths (`./inp/*.inp`) are relative to the CWD:

```bash
./main [jobname]              # legacy_final_vesion: reads legacy.in -> legacy.out (or jobname.out)
./w_asym [jobname]            # reads jobname.in -> jobname.out (default w_asym)
./w_pert < w_pert.in > out    # reads stdin, writes stdout
./get_yk_new <w_pert.out> <w_asym.out> <legacy_y.out> <R_Ai.txt> <out> <ZU|ZD|W+|W-> <NLO|NNLO> <pdf>
```

There is no unit-test suite. Validation is by `diff` against checked-in reference outputs: `w_pert_08112022/test_runs/`, `w_asym_08112022/test_runs/`, and `get_yk_new/grids_test_ykR/` + `grids_test_wp_tev2/` (see `get_yk_new/CLAUDE.md`).

## Cross-project things to know

- **Grid lockstep**: `w_pert`, `w_asym`, and Legacy must all be run on identical `(Q, qT, y)` grids (`inp/q_grid*.inp`, `qt_grid.inp`, `y_grid.inp`, same active-range indices). `get_yk_new` reads all three files row-by-row and aborts if any triplet differs by more than `1e-8`. When regenerating one input, regenerate the others with the same grid.
- **Legacy output format depends on `LTO`** (2nd field of line 2 of `legacy.in`; `-1`=LO only, `0`=full CSS+Y, `1`=ΔΣ integrated, `3`=Y-piece only): `LTO=3` is what `get_yk_new` expects; `LTO=0` writes an incompatible 4-column format.
- **What Legacy's Makefile actually builds**: `main.for`, `pert.for`, `pert_vj.for`, `res.for`, `pda.for`, `pion.for`, `CT14Pdf.f`, the five `*Pac02b.for`, `lhapdf.cpp`, linked with HOPPET + LHAPDF. `main_hj.for`, `pert_hj.for`, and `LhaPac02b.for` (link line commented out) are not built. `setup.sh` is a stale ATLAS `lsetup`/`/home/yfu` script, not the MSU HPCC environment — follow `README.md` instead.
- **Shared source copies**: `EvlPac02b/PrzPac02b/QcdPac02b/UtlPac02b/EwkPac02b.for`, `CT14Pdf.f`, `lhapdf.cpp` are duplicated, not shared, across Legacy, `w_pert`, and `w_asym` (`gauss.f`, `util.f` exist only in the latter two; Legacy uses `pda.for`/`pion.for` instead). A fix in one copy does not propagate. Current state (checked by md5): all are byte-identical across the three except two, where Legacy differs from `w_pert`/`w_asym` (which match each other):
  - `QcdPac02b.for` `FUNCTION NFL(AMU)`: Legacy picks nf=3/4/5 from `getthreshold` quark masses; `w_pert`/`w_asym` hardcode `nfl = 5` (old code left commented out). Intentional — don't "sync" it, but be aware fixed-order and resummed pieces use different nf schemes.
  - `UtlPac02b.for`: Legacy has two extra comment lines (`C_yfu change the MAXINT`); no code difference.
- **Physics-code lineage**: Legacy `main.for`/`res.for`/`pert.for` carry long dated change-log comments (`CCPY` = the original author's edits). `keep_codes/` in Legacy and `get_yk_new` holds timestamped snapshots/`.orig` files, and `*.f.03102021`, `asym_wrong_*.f`, `pert.for.delsig`, `pert.for.asympto` are historical variants — not built by any Makefile. Edit the file the Makefile names, not the variants.
- **Precompiled artifacts are checked in** (`get_yk_new_origin`, `get_yk_new_NoRAI`, `main`, `*.o`, `*.co`, `hoppet` binaries, `.pds`/`.LHgrid` PDF tables, `fort.*`). Don't edit or regenerate them by hand; `.gitignore` in Legacy excludes `main`, `*.o`, `fort.*`, `*.out`.
- **Parallelism**: Legacy is effectively **serial** (no OpenMP/MPI in the source, despite the README/SLURM templates setting `OMP_NUM_THREADS`). Long grids are parallelized by sharding along the outermost (Q) loop into a SLURM array: `make_shards.py` → `submit_*_array.sh` (`run_*_array.sb`) → `merge_shards.py`. Shard by Q only; splitting on y/qT interleaves rows and breaks the concatenation merge. The shard tooling is **not** in `legacy_final_vesion/`: it lives in `w_asym_08112022/` (for `w_asym`) and in `get_yk_new/E211/legacy/` (for Legacy); `get_yk_new/E211/resbos/` only has the array submit/run scripts (no make/merge). `submit_w_asym_array.sh JOBNAME 0-N[%throttle]` requires the range to start at 0.
- **`get_yk_new/E211/`** is an HPCC run snapshot (copies of each executable, `E211_*.in`, SLURM logs, shard directories) for one production campaign — treat as data/reference, not source. The top-level `E211_*` inputs in the sibling repos are the same campaign's steering files.
- Fixed-format Fortran: statements start in column 7, continuation in column 6, 72-column limit. Preserve this when editing; don't reflow.

## One-file campaign driver

`run_process.py` + `7TeV_WpWm_example.ini` (+ `shrds_scripts/`) prepare and chain a whole run from one INI: `python3 run_process.py 7TeV_WpWm_example.ini [--submit] [--reuse]`. The INI names a `source` folder (holding `get_yk_new/ legacy/ resbos/ w_asym/ w_pert/`; `templates/7TeV_WpWm/` is the reference one) and a new `dest`. It copies only what a run needs (executables, the `inp/*.inp` grids, `make_dummy_rai.py`; ~9 MB), rewrites the `.in` files, and writes `dest/submit_all.sh`, which chains SLURM jobs with `--dependency=afterok` by parsing the `Submitted merge job N` line printed by the shard submit scripts.

- **The source needs no outputs.** `get_yk_new` needs the `w_pert`/`w_asym`/Legacy(LTO=3) outputs and `resbos` needs the Legacy(LTO=0) and Yk grids, so each of those jobs stages its inputs (relative symlinks by default, real copies with `[paths] stage_outputs = copy`) when it *starts*, and stops with a clear error if one is missing. Nothing is linked or copied at prepare time; `Resbos_grids/` is created empty. The driver generates the `get_yk_new` and `resbos` job scripts itself (the source's `run_resbos.sb` is not used).
- **Executables** are taken from `bin/<folder>/<exe>` next to `run_process.py` by default — never from the template, whose binaries are stale snapshots. Override with `[executables]` (`dir = <dir>` with the `<dir>/<folder>/<exe>` layout, or one key per code: `w_pert`, `w_asym`, `legacy`, `get_yk_new`, `resbos`; a key wins over `dir`; `dir` may point at a template to opt in to its binaries). All five are checked before anything is copied, and a missing one aborts with a hint. `setup_resbos_legacy.sb` (plain bash, no python/conda; `sbatch` it or `bash -l` it on the HPCC, from this folder) copies the top-level sources of the five source dirs to `build/<folder>/`, does a from-scratch `make` with module/LHAPDF/HOPPET/ROOT settings taken from env vars (`MODULES`, `LHAPDF`, `HOPPET`, `ROOT_SETUP`; defaults = the HPCC paths of `[environment]`), and writes `bin/<folder>/<exe>` (`OUT=<dir>` changes the destination; `ONLY="legacy resbos"` builds a subset). It overrides Legacy's hardcoded `HOPPETLIBS` with `$HOPPET`. The source repos and their checked-in binaries are never touched. Tested locally only with a fake `gfortran` (no real compiler here).
- **Templates are the source's own W+ `.in` files** (named under `[templates]`); ECM, `lha_<pdf>`, `Type_V`/`JWTYPE` and the active range are patched by `> comment` marker, not line number. `[resbos]`/`[cuts.<run>]` keys and the legacy `bmax`/`nonpert` are optional overrides; with no `runs` key there is one run, `default`, that keeps the template's own cuts. W+ output is byte-identical to the Example's when the same template is used.
- **`[resbos] runs`** does not affect `w_pert`/`w_asym`/legacy/`get_yk_new` at all, so listing several names shares one upstream pass across all of them (one `get_yk_new`/Legacy-main dependency, several `resbos` jobs). `--resbos-only` (needs an already-prepared `dest` whose `legacy/*_main.out` and `get_yk_new/*.out` already exist) skips straight to (re)writing and, with `--submit`, submitting just the `resbos/*.in` + `run_*.sb` of `[resbos] runs` against those existing outputs — no dependency wait, and none of `w_pert`/`w_asym`/legacy/`get_yk_new` is touched. It re-copies the current `resbos_root` from `[executables]`, so a rebuilt `resbos` is picked up without a fresh full run.
- **Sharding uses `shrds_scripts/`** (the `w_asym` set). Per folder the driver rewrites the shard scripts with `w_asym` renamed to the code (`w_pert`, `w_asym`, `legacy`; Legacy's executable is `main`) and `HEADER_LINES` set to 3/3/15, and appends a row-count check to the merge job (Legacy `*_Y` = 3 lines/point, others 1). The shard scripts already present in the source are ignored. **`templates/7TeV_WpWm/legacy/merge_shards.py` has `HEADER_LINES = 3` (Legacy needs 15): merging with it silently produces a wrong file**; its other Legacy shard scripts are still the `w_asym` text (`w_asym_shard` job names/logs). Shard counts come from `[shards]` (1 is valid). **`resbos` is never sharded** (VEGAS MC).
- Legacy runs twice per boson (`_Y` = LTO 3 → `get_yk_new`; `_main` = LTO 0 → `resbos`); `resbos` runs once per entry of `[resbos] runs`.
- Guards: the three grid sets (`w_pert`/`w_asym`/`legacy` `inp/`) must be md5-identical; job names >95 chars are refused (Fortran `character*100`); `--reuse` refuses to reuse a `<job>_shards/` dir whose `.in` changed (`submit_*_array.sh` would silently keep the old shards, i.e. the old energy).
- Only `W+`/`W-` are wired up (the Z chain through `get_yk_new`/`resbos` is unconfirmed). `get_yk_new` at `NLO` gets an all-ones dummy `R_Ai`; `NNLO` needs `[get_yk_new] r_ai_<tag>`.
- Tested locally: input generation from both `Example_E211_ATLAS8TeV_WpWm` and `templates/7TeV_WpWm`, script syntax, `submit_all.sh` against a fake `sbatch` (ids/dependencies), and the merge/`get_yk_new`/`resbos` job scripts run end to end on a 12-point grid with fake executables. The real SLURM/HPCC run has not been done (no `sbatch`/`gfortran` here). `get_yk_new` walltime (12 h) and memory (16 G) are guesses.
