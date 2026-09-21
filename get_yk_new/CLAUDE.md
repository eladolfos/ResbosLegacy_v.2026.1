# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Fortran 77 legacy physics code that computes K-factors for QCD calculations. The main program `get_yk_new` combines perturbative and asymptotic results from multiple sources to produce refined kinematic K-factors for particle production processes (DY, W±, Z boson production).

**Physics Context**: The code processes output from:
- Perturbative QCD calculations (from MCFM)
- Asymptotic soft-collinear resummation results
- Angular function corrections (R_Ai from FCFM)
- Legacy baseline calculations from ResBoS

Output K-factors (y and ykR) are used in particle physics phenomenology for accurate cross-section predictions.

## Building and Running

### Compile
```bash
make
# Produces executable: ./get_yk_new
```

### Run (ZU particles with NLO order, LHC 8 TeV)
```bash
./get_yk_new w_pert_ZU.out w_asym_ZU.out legacy_y_ZU.out R_Ai_LHC8_66Q116.txt output_ykR_ZU.out.NLO ZU NLO CT14nn.00
```

### Run (W+ production with NNLO order, Tevatron)
```bash
./get_yk_new w_pert_wp_tev2.out w_asym_wp_tev2.out legacy_y_wp_tev2.out R_Ai_tev2_wp.txt output_ykR_wp_tev2.out.NNLO W+ NNLO CT14nn.00
```

### Clean
```bash
make clean
```

## Code Architecture

### Main Program: `get_yk_new.f` (~1564 lines)

**Core Function**: Reads 4-5 input files and produces K-factor output by:
1. Reading perturbative results (w_pert file) with L0 and A3 terms
2. Reading asymptotic/resummed results (w_asym file) 
3. Reading legacy baseline calculations (legacy_y file)
4. Applying angular function corrections (R_Ai values)
5. Combining results according to QCD order (NLO vs NNLO)

**Key Parameters**:
- `Boson`/`itype`: particle type, checked verbatim against `'ZU'`, `'ZD'`, `'W+'`, `'W-'` (get_yk_new.f:181-185) — other strings (e.g. `DY`) fall through without a matching branch
- `iorder`: calculation order (NLO or NNLO)
- `iYPiece`: controls which calculation piece to use
- `iYPSw`: perturbative calculation switch
- `iKFacP`, `iKFacY`, `iYGrid`: K-factor and grid switches

**Upstream Input Generation**: `w_pert_*.out`, `w_asym_*.out`, and `R_Ai_*.txt` are not produced by this repo — they come from separate MCFM, ResBos2, and FCFM runs (see the `.f` variants under `grids_test_ykR/compare_old_pert_asym_out/codes_modified/` and `grids_test_wp_tev2/codes_revised/` for the perturbative/asymptotic generator code actually used to build the checked-in test grids, and the `inp/q_grid*.inp`, `inp/qt_grid.inp`, `inp/y_grid.inp` files for the grid definitions).

**Input File Format** (tab/space-separated numerical grids):
- `w_pert_*.out`: Header with metadata, then Q, qT, y columns with perturbative terms (pertL0_1, pertA3_1, sums)
- `w_asym_*.out`: Similar structure with asymptotic terms (asymL0_1, asymA3_1, sums)
- `legacy_y_*.out`: Header with physics parameters, then baseline calculations
- `R_Ai_*.txt`: Header with ECM and type, then Q, qT, y with angular correction factors (R_A0 through R_A3)

### Directory Structure

- **Root**: Main source code and standard test data files (also used as the default working directory for `make`/`./get_yk_new` runs)
- **check_R_Ai/**: Verification test case for angular function corrections (uses CT14.00 PDF); contains its own copy of `get_yk_new.f` plus `.in`/`.out` reference pairs
- **grids_test_ykR/**: Test directory for ZU particles (LHC 8 TeV) - contains both input grids and reference output, plus `compare_old_pert_asym_out/codes_modified/` (older `asym.f`/`pert.f` variants used to cross-check upstream MCFM/ResBos output)
- **grids_test_wp_tev2/**: Test directory for W+ production (Tevatron) - alternate grid for validation, plus `codes_revised/` (revised `asym.f`/`get_yk_new.f`/`w_pert.f` used to regenerate that grid)
- **keep_codes/**: Historical, timestamped versions of `get_yk_new.f` (backups from Oct 2020 - Mar 2022, e.g. `get_yk_new.f.03182022`)

### Key Source Variants

- **get_yk_new.f**: Current/main version (as of Mar 2022) — this is the only source file built by the `Makefile`
- **get_yk_new_NoRAI**, **get_yk_new_origin**: Precompiled ELF binaries checked into the repo (not source — `file <name>` confirms), kept as reference executables from before the R_Ai angular-correction update. Do not edit; if you need their behavior, check out the corresponding revision from `keep_codes/` and rebuild instead.
- **keep_codes/**: Timestamped source snapshots for reference/diffing against the current `get_yk_new.f`

## Compilation Flags

**Current**: `-O3 -fno-automatic`
- `-O3`: Full optimization for production code
- `-fno-automatic`: Force static allocation (legacy Fortran compatibility)

**Alternative (commented)**: `-g -fvxt` for debugging with extra checks

## Important Notes for Development

1. **Fortran 77 Legacy Code**: Uses fixed-format Fortran with column-based syntax. Preserve spacing and formatting when editing.

2. **Input File Dependencies**: The program requires all 4-5 input files to be present and properly formatted. Changes to input parsing must be verified against all test cases.

3. **Test Data**: Use test directories (`grids_test_ykR/`, `grids_test_wp_tev2/`) to validate changes. Legacy output files show expected results.

4. **Physics Order Handling**: NLO and NNLO calculations use different combination schemes for perturbative and asymptotic terms. The `iorder` parameter controls which branches are executed.

5. **Grid Alignment**: K-factor values are grid-aligned with input (Q, qT, y triplets). Ensure alignment is preserved when modifying output logic.

6. **PDF Set Parameter**: Last argument (e.g., `CT14nn.00`) is typically a PDF set identifier - may be used for logging or validation.

## Common Tasks

**Verify a change produces correct output**:
```bash
make clean && make
./get_yk_new w_pert_ZU.out w_asym_ZU.out legacy_y_ZU.out R_Ai_LHC8_66Q116.txt test_output.out ZU NLO CT14nn.00
# Compare test_output.out with legacy_ykR_ZU.out.NLO
diff test_output.out legacy_ykR_ZU.out.NLO
```

**Test both physics processes**:
```bash
# Test ZU (with both NLO and NNLO)
./get_yk_new w_pert_ZU.out w_asym_ZU.out legacy_y_ZU.out R_Ai_LHC8_66Q116.txt test_ZU_NLO.out ZU NLO CT14nn.00
./get_yk_new w_pert_ZU.out w_asym_ZU.out legacy_y_ZU.out R_Ai_LHC8_66Q116.txt test_ZU_NNLO.out ZU NNLO CT14nn.00

# Test W+ (Tevatron)
./get_yk_new w_pert_wp_tev2.out w_asym_wp_tev2.out legacy_y_wp_tev2.out R_Ai_tev2_wp.txt test_WP.out W+ NLO CT14nn.00
```
