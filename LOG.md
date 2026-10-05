# LOG

Investigation and code history for the resbos-legacy kinematics and resNLO work (Oct 2026).
Records what was found, with the evidence, what was **refuted**, and what changed in the code.
Commit hashes are on `main`.

---

## 1. Summary

Two separate failures were mixed up during the investigation. They have different causes.

| # | Failure | Where | Cause | Status |
|---|---|---|---|---|
| A | Legacy `ERROR IN PERTURB ... Stopping in QUIT` kills a shard | `legacy_final_vesion/pert.for` | Grid point with Q > ECM (parton x > 1) | Handled: active Q capped at Q ≤ ECM (`run_process.py` check) |
| B | resNLO (and nloasy) `TOTAL CROSS SECTION = NaN`, `exit code 0` | `resbos/resbosm.f`, `YMAXIMUM` | Sampled (Q, qT) cells with no phase space: `acosh(X)` with X < 1 | Fixed by a guard (`6f44c65`), validated at 1800 GeV and at 6 energies |

Failure B was first blamed on `Q == ECM` exactly. That hypothesis was **refuted** (see §4).

---

## 2. Findings

### 2.1 Legacy PERTURB (failure A)

`legacy_final_vesion/pert.for`, subroutine around lines 466–526:

```
TM_V  = sqrt(QT_V^2 + Q_V^2)
RTAUP = (TM_V + QT_V) / ECM
X1LOW = RTAUP * exp(Y_V)
X2LOW = RTAUP / exp(Y_V)
```

- If `X1LOW > 1` and `X2LOW > 1`, legacy writes `ERROR IN PERTURB` and calls `QUIT`, which kills the whole shard.
- At `qT → 0`, `y = 0` this is `Q > ECM`.
- Observed on ppbar, ECM = 1800 GeV: every Q-shard with Q ≤ 1800 finished clean; every Q-shard with Q ≥ 1900 died within seconds (`x1Low = 1.0556`, which is 1900/1800).
- `Q == ECM` passes legacy (`RTAUP = 1.0` is not `> 1.0`, and shard 11 with Q = 1800 ran clean).

### 2.2 Legacy output does not depend on the "Full range" line

- `main.for` reads only the `Active setting` line (`IQTMN, ..., IQST`, around line 715), then closes the input file. The `Full range` line is never read.
- Changing `Full range` cannot change any value in the grid output.

### 2.3 resbos sizes its grids from the data, not from a header

- `PreReadIn` / `Sample` (`resbos_root.f:4380–4525`) count `iD1, iD2, iD3` and `iY1, iY2, iY3` from the file itself.
- A truncated grid (fewer Q points) is sized correctly. Array bounds are not the cause.

### 2.4 The resbos NaN (failure B)

`resbos/resbosm.f`, `YMAXIMUM`, line ~5597:

```fortran
ACOSH(X) = DLOG(X + DSQRT(X**2 - 1.0))
X = (ECM^2 + Q^2) / (2 * ECM * sqrt(Q^2 + qT^2))   ! = cosh(y_max)
YMAX = ACOSH(X)                                     ! NaN when X < 1
```

- `WRESPH` (used by `NPART = 11` resummed and `NPART = 12` nloasy) calls `YMAXIMUM`, then computes `YBOOST`, `X1 = sqrt(tau)*exp(YBOOST)`, `X2 = tau/X1`, and the weight. One NaN sample makes the VEGAS sums NaN from iteration 1.
- `WLOPH` (`NPART = 13`, nlodsi) also calls `YMAXIMUM` and did not produce NaN in the same windows.
- `hvvresph.f`, `AARESPH` and `WGRIDPH` also call `YMAXIMUM`. They are **not validated**.

Exact kinematic boundary, derived from momentum conservation for a boson (mT, qT, y) plus one recoil parton with transverse momentum qT:

```
x1,2 = (mT·e^±y + qT·e^±yk) / √S,   mT² = Q² + qT²,   S = ECM²
existence  ⇔  cosh(y) ≤ (S + Q²) / (2·√S·mT)
at y = 0:   qT ≤ (S − Q²) / (2·√S)
```

Forbidden (Q, qT) cells, counted with this criterion over the full box:

| Energy | Window | Forbidden cells | Observed |
|---|---|---|---|
| 1800 GeV | Q ≤ 1700, full qT | 2084 / 7191 | NaN (all three resNLO and nloasy) |
| 1800 GeV | Q ≤ 1000, qT ≤ 600 | 0 / 4280 | finite, 1251.349 pb |
| 7 TeV | Q ≤ 1700, full qT | 0 / 7191 | finite |
| 7 TeV | Q ≤ 5000, full qT | 251 / 12240 | finite (not explained, see §5) |
| 13 TeV | Q ≤ 5000, full qT | 0 / 12240 | finite |
| 5 TeV | Q ≤ 4900, full qT | 1541 / 12087 | 268 KB stub, no log (see §5) |

Legacy writes zeros in those cells, so the grid files contain no NaN. The NaN is created inside resbos.

### 2.5 Input grids are clean

The NaN is not present in the inputs. For Q ≤ 1700 at 1800 GeV: `legacy_main.out` = 1028328 lines = 15 + 47·143·153, and `Yk.out` = 3084956 lines = 17 + 3·47·143·153. No `NaN` or `Inf` tokens.

### 2.6 Validated results after the guard

All runs are W+ unless noted. Relative difference is `(resNLO − (nloasy + nlodsi)) / (nloasy + nlodsi)`.

| Energy | Window | resNLO (pb) | nloasy + nlodsi (pb) | Rel. diff |
|---|---|---|---|---|
| 1800 GeV | Q ≤ 1000, qT ≤ 600 | 1251.349 ± 0.21 (regression, resNLO only) | — | — |
| 1800 GeV | Q ≤ 1700, full qT | 1251.318 ± 0.22 | 1149.649 + 99.002 = 1248.65 | +0.21 % |
| 1960 GeV | Q ≤ 1900, full qT | 1368.686 ± 0.24 | 1289.221 + 78.951 = 1368.17 | +0.04 % |
| 5 TeV | Q ≤ 4900 | 4228.02 ± 1.06 | 5041.29 − 746.47 = 4294.82 | −1.56 % |
| 7 TeV | full | 6029.50 ± 1.48 | 7684.90 − 1538.92 = 6145.99 | −1.90 % |
| 8 TeV | full | 6920.27 ± 1.69 | 9037.63 − 1976.03 = 7061.60 | −2.00 % |
| 13 TeV | full | 11273.69 ± 2.68 | 15957.39 − 4413.88 = 11543.51 | −2.34 % |

The regression run (Q ≤ 1000, qT ≤ 600) reproduces 1251.349 pb exactly with the guard, so the guard does not change the result in a window without forbidden cells.

Note: the Q ≤ 1700 run did **not** come out larger than the Q ≤ 1000 run (1251.318 vs 1251.349). The added region (Q in (1000, 1700], qT > 600) contributes less than the MC error. That is consistent, not a bug.

Before the fix, the `NLO.root` (`hadd` of nlodsi and nloasy) looked fine (1 GB) only because nloasy was a stub. Its real content was only the nlodsi piece (99 pb). Any earlier NLO total from that campaign is invalid.

---

## 3. Code changes

| Commit | File(s) | What it does |
|---|---|---|
| `6f44c65` | `resbos/resbosm.f` | `YMAXIMUM`: if `X < 1`, set `YMAX = 0` (zero phase space, zero weight) instead of `acosh(X)` = NaN. Points with `X ≥ 1` are unchanged. |
| `2ae418e` | `run_process.py` | `check_q_kinematics` and `check_q_kinematics_points` reject active Q above ECM (the PERTURB limit, `Q ≤ ECM`). Messages and docstrings describe only that failure. |
| `02b72d9` | `Examples/ppbar_WpWm_CT25_Tevatron_1800*`, `…_1960…`, `…_1800_1960…`, `…_only_Q1000`, `…_regression_guard`, `…_Q1700_fullqT_guard` | Tevatron Run I / Run II examples with the validated caps (47 and 49), plus the diagnostic, regression and physics-change tests. |
| `db605d6`, superseded by `5d65638` | `scripts/validate_campaign.py`, `run_process.py` | First version of the post-run check, run from a background `finalize_campaign.sh`. |
| `5d65638` | `run_process.py`, `scripts/validate_campaign.py` | `write_analysis()` writes `analysis.sb`, a single-core SLURM job submitted last by `submit_all.sh`. It waits until no pending/running job with the campaign's name prefix remains (ignoring `DependencyNeverSatisfied`), runs `validate_campaign.py`, and appends TOTAL TIME and the tables to `timing.log`. The old background finalize was removed: it died with the login-node session. |
| `28e6642` | `Examples/validation_<E>_Wp_CT25_resNLO_NLO.ini` (1800, 1960, 5, 7, 8, 13 TeV) | Validation set, W+ only, `compute = resNLO, NLO`. |
| `f7c600d` | `Examples/production_<E>_WpWmZ0_CT25_resNLO_NLO.ini` | Production set, W+, W−, Z0, same windows, job names prefixed `prod`. |

### 3.1 `scripts/validate_campaign.py` (current behaviour)

Checks, per resbos piece (latest job per piece):

1. `TOTAL CROSS SECTION` is finite, the log says `finished with exit code 0`, and there are no `NaN` tokens.
2. Each `.root` in the dest is larger than `--min-root-mb` (default 100 MB). Smaller means a stub.
3. The merged grids (`legacy`, `w_pert`, `w_asym`, `get_yk_new`, `resbos/Resbos_grids`) contain no NaN/Inf. This scan is slow (~3 GB).
4. If resNLO, nloasy and nlodsi all exist: `|resNLO − (nloasy+nlodsi)| / (nloasy+nlodsi)` above `--cross-tol` (default 3 %) is a WARN, not an ERROR.

Output: two tables (pieces; resNLO vs sum) and `VALIDATION PASSED|FAILED`, written to `validation_summary.txt` and printed.

Verified: PASS on the Q ≤ 1000 and Q ≤ 1700 runs; FAIL with 24 errors on the old NaN campaign (NaN cross sections and 0.06 MB stubs).

### 3.2 Constraints the `.ini` files now follow

- `active` Q upper index ≤ ECM. Tevatron 1800 → index 47 (1700 GeV). Tevatron 1960 → index 49 (1900 GeV). 5 TeV → index 79 (4900 GeV, one step below the grid's last point, which equals ECM). 7, 8, 13 TeV → full grid (index 80, 5000 GeV).
- The 13 TeV template's y grid has 155 points, not 143: `active` y must be `1 155 1`.
- The `Tevatron_WpWm` template is p-pbar: `ibeam = -1`. pp templates use `ibeam = 1`.

---

## 4. Refuted or retracted claims

Do not reuse these.

| Claim | Status | Why |
|---|---|---|
| "NaN appears when active Q_max equals ECM exactly" | **Refuted** | With Q_max = 1700 (100 GeV margin) the NaN persisted (run `18358263`, `Maximal invariant mass (Q): 1700`). |
| "Strict `Q < ECM` is required" (`run_process.py` at `2ae418e`'s predecessor) | **Retracted** | Based on the refuted hypothesis. Restored to `Q ≤ ECM` in `2ae418e`. |
| "The stencil in `MatchGrids1/2` reads past the truncated grid" | **Refuted** | `PreReadIn`/`Sample` size the grids from the data (§2.3). |
| "The `Full range` mismatch corrupts legacy output" | **Refuted** | Legacy never reads that line (§2.2). |
| "NLO works, only resNLO fails" | **Retracted** | nloasy was also NaN. `NLO.root` was a `hadd` with a stub (§2.6). |
| "The 7 TeV full grid agrees with the forbidden-cell rule" | **Open** | 251 forbidden cells at Q ≤ 5000 and the run was finite. See §5. |

---

## 5. Open questions

1. **7 TeV full grid**: 251 forbidden cells (Q ≤ 5000, qT > ~1700) and a finite result. Either those cells are never sampled, or they are zero-weighted some other way. Not verified.
2. **5 TeV stub** (pp Z0 campaign, `pp_Z0_CT25_5_7_8_13TeV_resNLO_NLO`): stub `resNLO.root`, no log. Consistent with forbidden cells (1541 at Q ≤ 4900), but the job may never have run. Not reproduced with the guard. The 5 TeV validation run (`validation_5TeV`) gives finite resbos results for all three pieces, and its resbos-level checks passed. Its full grid NaN scan had not finished when this was written, so the validation verdict is not yet confirmed.
3. **Energy trend**: the resNLO vs NLO difference goes from +0.04 % (1960) to −2.3 % (13 TeV), with a sign change between 1960 and 5 TeV. Possible causes: the qT ≤ 3000 GeV and Q ≤ 5000 GeV cut in the resummed grid versus the fixed-order pieces, or a scale or PDF mismatch. Not investigated.
4. **Other samplers**: `hvvresph.f`, `AARESPH`, `WGRIDPH` also call `YMAXIMUM` and are covered by the guard in code, but none of them has been run with it.
5. **W− and Z0 validation**: the validation set is W+ only. Production runs W+, W−, Z0 with the same windows, which have not been validated individually.
6. **Guard semantics**: a zero-weight event still appears in the ntuple. Cross sections are unaffected. Unweighted event files would contain zero-weight rows.

---

## 6. How to run a campaign with these fixes

1. Recompile resbos with the guard: `sbatch setup_resbos_legacy.sb` (ONLY="resbos" is enough if the others are built).
2. `git pull` on the HPCC so `run_process.py` and `scripts/validate_campaign.py` are current.
3. Prepare and submit: `python3 run_process.py Examples/<file>.ini --submit`.
4. When the queue empties, `analysis.sb` runs the validator. Read `timing.log` (tables at the end) or `validation_summary.txt`.
5. Accept a campaign only if `VALIDATION PASSED` and every resNLO vs (nloasy+nlodsi) difference is within the expected trend (§5.3).

For campaigns already run before `5d65638`, run the validator by hand:

```bash
sbatch --partition=general-long --time=02:00:00 --mem-per-cpu=4G \
  --wrap="python3 scripts/validate_campaign.py <dest>"
```

---

## 7. Physics reference

- **CSS resummation**: `dσ/dqT` at small qT is the Fourier-Bessel integral over b of the Sudakov-resummed `W(b, Q, x1, x2)` (`e^{-S(Q,b)}` times the PDFs and coefficients). It is valid at small qT only.
- **Y-term matching**: `resNLO = CSS + Y`, where `Y = (fixed order) − (asymptotic)`. The matching reproduces fixed-order NLO at all qT, so `resNLO ≈ nloasy + nlodsi` in total cross section up to O(α_s²) and MC error. That is the consistency check in §3.1.
- **Kinematic limits**: §2.4. Any sampled `(Q, qT, y)` with `cosh(y) > (S+Q²)/(2√S·mT)` has no phase space, so its contribution is zero by PDF support (`f(x) = 0` for `x > 1`).
