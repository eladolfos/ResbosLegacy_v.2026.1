# get_yk_new: K-Factor Calculator for QCD Processes

A Fortran 77 legacy physics code that computes kinematic K-factors for Drell-Yan (DY), W±, and Z boson production by combining perturbative QCD calculations, asymptotic resummation results, and angular correction functions.

## Table of Contents

- [Overview](#overview)
- [Prerequisites](#prerequisites)
- [HPCC Installation & Setup](#hpcc-installation--setup)
- [Building the Code](#building-the-code)
- [Using the Program](#using-the-program)
- [Process-Specific Workflows](#process-specific-workflows)
- [Grid consistency between w_pert, w_asym, and legacy_y](#grid-consistency-between-w_pert-w_asym-and-legacy_y)
- [Running without a real R_Ai file](#running-without-a-real-r_ai-file)
- [File Organization](#file-organization)
- [Troubleshooting](#troubleshooting)

## Overview

The `get_yk_new` program processes output from multiple QCD calculation codes (MCFM, ResBos, FCFM) to compute refined K-factors. It combines:
- **Perturbative results** (leading order L0, next-to-leading A3 terms)
- **Asymptotic/resummation results** (soft-collinear resummation)
- **Angular correction functions** (R_Ai factors from FCFM)
- **Baseline calculations** (legacy results)

The K-factors are essential for phenomenological predictions of particle production cross sections at hadron colliders.

## Prerequisites

### Required Software
- **Fortran 77 Compiler**: `gfortran` (included in most HPC environments)
- **LHAPDF**: PDF library (for PDF handling)
- **ROOT**: For analysis and data handling (optional but recommended)
- **MCFM**: QCD calculation code (for generating perturbative input)
- **ResBos2**: Resummation code (for asymptotic results)
- **HOPPET**: PDF evolution code (optional)

### Knowledge
- Basic familiarity with QCD calculations
- Understanding of K-factors in particle phenomenology
- HPC cluster usage (SLURM job submission, module system)

## HPCC Installation & Setup

### Step 1: Set Up Your Environment

Add the following to your `~/.bashrc` file on HPCC:

```bash
# ============================================
# get_yk_new environment setup (HPCC/MSU)
# ============================================

# Test for profile.dos compatibility
test -f /etc/profile.dos && . /etc/profile.dos

# LHAPDF configuration
export LD_LIBRARY_PATH=/mnt/home/lopezels/InstallSources/LHAPDF/lib:$LD_LIBRARY_PATH
export PATH=/mnt/home/lopezels/InstallSources/LHAPDF/bin:$PATH

# HOPPET configuration
export LD_LIBRARY_PATH=/mnt/home/lopezels/InstallSources/HOPPET1/lib:$LD_LIBRARY_PATH
export PATH=/mnt/home/lopezels/InstallSources/HOPPET1/bin:$PATH

# MCFM path (perturbative calculations)
export PATH=/mnt/home/lopezels/InstallSources/MCFM-10.3:$PATH

# ROOT configuration (must come before ResBos2)
source /mnt/home/lopezels/InstallSources/ROOT/bin/thisroot.sh

# ResBos2 configuration (resummation code)
export PATH=/mnt/home/lopezels/InstallSources/ResBos2/bin:$PATH
export LD_LIBRARY_PATH=/mnt/home/lopezels/InstallSources/ResBos2/lib:$LD_LIBRARY_PATH

# OpenMP settings for parallel processing
export OMP_NUM_THREADS=28        # Adjust based on allocation
export OMP_STACKSIZE=2G          # Stack size for OpenMP threads

# Load aliases if they exist
test -s ~/.alias && . ~/.alias
```

### Step 2: Load HPCC Modules

On HPCC, load required modules:

```bash
module purge
module load GCC/11.2.0  # or your preferred compiler module
module load LHAPDF      # if available as a module
```

Check available modules:
```bash
module avail
```

### Step 3: Clone or Copy the Repository

```bash
# Option 1: Copy from existing location
cp -r /path/to/get_yk_new ~/MyResearch/get_yk_new
cd ~/MyResearch/get_yk_new

# Option 2: Initialize in your working directory
mkdir -p ~/MyResearch/get_yk_new
cd ~/MyResearch/get_yk_new
# Copy files (see below)
```

### Step 4: Verify Installation

```bash
# Source your updated bashrc
source ~/.bashrc

# Check Fortran compiler
gfortran --version

# Compile the code
make clean
make

# Verify executable was created
ls -lh get_yk_new
```

## Building the Code

### Quick Compile
```bash
cd ~/MyResearch/get_yk_new
make
```

### With Debug Symbols
Edit `Makefile` to use debug flags:
```makefile
COMPFLAGS = -g -fvxt -fno-automatic
LINKFLAGS = -fvxt -g -fno-automatic
```
Then compile:
```bash
make clean
make
```

### Clean Build
```bash
make clean
make
```

## Using the Program

### Basic Syntax

```bash
./get_yk_new <w_pert_file> <w_asym_file> <legacy_y_file> <R_Ai_file> <output_file> <boson_type> <order> <pdf_set>
```

### Parameter Description

| Parameter | Description | Examples |
|-----------|-------------|----------|
| `w_pert_file` | Perturbative QCD results from MCFM | `w_pert_ZU.out`, `w_pert_wp_tev2.out` |
| `w_asym_file` | Asymptotic/resummation results | `w_asym_ZU.out`, `w_asym_wp_tev2.out` |
| `legacy_y_file` | Baseline calculation (ResBoS output) | `legacy_y_ZU.out`, `legacy_y_wp_tev2.out` |
| `R_Ai_file` | Angular correction functions | `R_Ai_LHC8_66Q116.txt`, `R_Ai_tev2_wp.txt` |
| `output_file` | Output K-factor file | `legacy_ykR_ZU.out.NLO` |
| `boson_type` | Particle type | `ZU`, `ZD`, `W+`, `W-`, `DY` |
| `order` | Calculation order | `NLO`, `NNLO` |
| `pdf_set` | PDF set identifier | `CT14nn.00`, `CT18nnlo.00` |

### Example Runs

#### Example 1: ZU Drell-Yan at NLO (LHC 8 TeV)
```bash
./get_yk_new w_pert_ZU.out w_asym_ZU.out legacy_y_ZU.out R_Ai_LHC8_66Q116.txt output_ykR_ZU.out.NLO ZU NLO CT14nn.00
```

#### Example 2: ZU Drell-Yan at NNLO
```bash
./get_yk_new w_pert_ZU.out w_asym_ZU.out legacy_y_ZU.out R_Ai_LHC8_66Q116.txt output_ykR_ZU.out.NNLO ZU NNLO CT14nn.00
```

#### Example 3: W+ Production at Tevatron (NLO)
```bash
./get_yk_new w_pert_wp_tev2.out w_asym_wp_tev2.out legacy_y_wp_tev2.out R_Ai_tev2_wp.txt output_ykR_wp_tev2.out.NLO W+ NLO CT14nn.00
```

#### Example 4: W+ Production at Tevatron (NNLO)
```bash
./get_yk_new w_pert_wp_tev2.out w_asym_wp_tev2.out legacy_y_wp_tev2.out R_Ai_tev2_wp.txt output_ykR_wp_tev2.out.NNLO W+ NNLO CT14nn.00
```

## Process-Specific Workflows

### Drell-Yan (DY) / Z Production

**For LHC 8 TeV Analysis:**
```bash
# Copy required files to your working directory
cp grids_test_ykR/w_pert_ZU.out .
cp grids_test_ykR/w_asym_ZU.out .
cp grids_test_ykR/legacy_y_ZU.out .
cp grids_test_ykR/R_Ai_LHC8_66Q116.txt .

# Run NLO calculation
./get_yk_new w_pert_ZU.out w_asym_ZU.out legacy_y_ZU.out R_Ai_LHC8_66Q116.txt my_ykR_ZU_NLO.out ZU NLO CT18nnlo.00

# Run NNLO calculation
./get_yk_new w_pert_ZU.out w_asym_ZU.out legacy_y_ZU.out R_Ai_LHC8_66Q116.txt my_ykR_ZU_NNLO.out ZU NNLO CT18nnlo.00

# Verify against reference
diff my_ykR_ZU_NLO.out grids_test_ykR/legacy_ykR_ZU.out.NLO
diff my_ykR_ZU_NNLO.out grids_test_ykR/legacy_ykR_ZU.out.NNLO
```

**For Custom Energy/Grids:**
You need to generate input files using MCFM, ResBos, and FCFM:
1. Run MCFM for perturbative calculations → produces `w_pert_*.out`
2. Run ResBos for asymptotic resummation → produces `w_asym_*.out`
3. Prepare baseline calculation → produces `legacy_y_*.out`
4. Generate R_Ai angular functions → produces `R_Ai_*.txt`

### W Boson Production

**For Tevatron (sqrt(s) = 1.96 TeV):**
```bash
# Copy required files
cp grids_test_wp_tev2/w_pert_wp_tev2.out .
cp grids_test_wp_tev2/w_asym_wp_tev2.out .
cp grids_test_wp_tev2/legacy_y_wp_tev2.out .
cp grids_test_wp_tev2/R_Ai_tev2_wp.txt .

# Run for W+ at NLO
./get_yk_new w_pert_wp_tev2.out w_asym_wp_tev2.out legacy_y_wp_tev2.out R_Ai_tev2_wp.txt ykR_wp_NLO.out W+ NLO CT14nn.00

# Run for W+ at NNLO
./get_yk_new w_pert_wp_tev2.out w_asym_wp_tev2.out legacy_y_wp_tev2.out R_Ai_tev2_wp.txt ykR_wp_NNLO.out W+ NNLO CT14nn.00

# For W-, just change boson type
./get_yk_new w_pert_wp_tev2.out w_asym_wp_tev2.out legacy_y_wp_tev2.out R_Ai_tev2_wp.txt ykR_wm_NLO.out W- NLO CT14nn.00
```

### Batch Processing on HPCC

**Create `run_kfactors.slurm`:**
```bash
#!/bin/bash
#SBATCH --job-name=kfactors
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --output=kfactors_%j.log

module purge
module load GCC/11.2.0

cd ~/MyResearch/get_yk_new

# Compile if needed
make

# Run multiple processes
echo "Computing ZU NLO K-factors..."
./get_yk_new w_pert_ZU.out w_asym_ZU.out legacy_y_ZU.out R_Ai_LHC8_66Q116.txt ykR_ZU_NLO.out ZU NLO CT18nnlo.00

echo "Computing ZU NNLO K-factors..."
./get_yk_new w_pert_ZU.out w_asym_ZU.out legacy_y_ZU.out R_Ai_LHC8_66Q116.txt ykR_ZU_NNLO.out ZU NNLO CT18nnlo.00

echo "Computing W+ NLO K-factors..."
./get_yk_new w_pert_wp_tev2.out w_asym_wp_tev2.out legacy_y_wp_tev2.out R_Ai_tev2_wp.txt ykR_wp_NLO.out W+ NLO CT14nn.00

echo "All K-factor calculations completed!"
```

**Submit the job:**
```bash
sbatch run_kfactors.slurm

# Check status
squeue -u $USER

# View output
tail -f kfactors_*.log
```

## Grid consistency between w_pert, w_asym, and legacy_y

`get_yk_new` reads `w_pert_*.out`, `w_asym_*.out`, and `legacy_y_*.out` in lockstep, one grid point per iteration, and requires their `(Q, qT, y)` triplet to match at every single row (`get_yk_new.f:170-179`, tolerance `1e-8`). The instant one file's row disagrees with the others, the program stops immediately:

```
Kinematical parameters do not match at line 16
x1p, x2p, x3p =   60.0000000       9.10599977E-02  -3.90000010
x1a, x2a, x3a =   80.0000000       9.10599977E-02  -3.90000010
X1R, X2R, X3R =   10.0000000       9.10999998E-02  -5.00000000
```

If you hit this, the three input files were generated over different grids and something upstream needs to be fixed — `get_yk_new` itself cannot reconcile mismatched grids.

### It's rarely the grid *files* that differ

`q_grid.inp`, `qt_grid.inp`, and `y_grid.inp` are typically copied verbatim into each generator's `inp/` directory (`w_pert/inp/`, `w_asym/inp/`, `legacy/inp/`, etc.) and stay byte-identical — `diff` them first to confirm/rule this out. The actual mismatch is almost always in the **`.in` config file** each generator (`w_pert`, `w_asym`, the legacy code) reads before running, which controls the grid two different ways:

1. **Which grid file is referenced.** Each `.in` file has a "grid file name" line per dimension, e.g.:
   ```
   ./inp/q_grid.inp                          > Q grid file name
   ./inp/qt_grid.inp                         > qT grid file name
   ./inp/y_grid.inp                          > y grid file name
   ```
   It's easy for one generator to point at a different file than the others — e.g. a reduced `q_grid_test.inp` (a handful of Q values, for quick test runs) instead of the full production `q_grid.inp`.

2. **Which subset of that file is actually computed.** Even when all three `.in` files reference the *same* grid file, an "Active setting" line selects which indices of it get computed, as `<start> <end> <step>` triplets for qT, y, and Q in that order:
   ```
   1  153  1   1  143  1  1 80  1            > Active setting
   ```
   This example means: qT indices 1–153 step 1 (all of them), y indices 1–143 step 1 (all of them), Q indices 1–80 step 1 (all of them) — i.e. the full grid. A narrower setting like `1 105 1  12 132 12  1 7 2` only computes qT indices 1–105, y indices 12,24,...,132, and Q indices 1,3,5,7 — a real subset, even though it's reading from the exact same grid file. Each `.in` file usually has its own "Full range" comment line documenting what the all-inclusive setting looks like for that file's grid — use it as the reference when re-aligning a mismatched generator, but double check it was actually updated when the grid file itself changed (it's just a comment, nothing enforces it stays in sync).

**To fix a mismatch**: diff the "grid file name" and "Active setting" lines across all the `.in` files feeding a single `get_yk_new` run, and make them identical (same referenced files, same start/end/step for qT, y, and Q), then regenerate whichever outputs used the wrong or narrower grid.

## Running without a real R_Ai file

You need MCFM/ResBos2/FCFM to produce a genuine `R_Ai_*.txt` file, but `get_yk_new` still requires an `R_Ai` argument even if you're not ready to generate one. **For NLO runs this is fine to fake** — the angular-correction values are never applied at NLO, so a placeholder file works and produces identical output to a real one. **For NNLO runs, faking it means you are explicitly turning off the angular correction**, not approximating it — only do this if that's actually what you want (see below).

### Why this works

In `get_yk_new.f`, the R_Ai correction factors are initialized to `1.0` and only overwritten by real interpolated values inside an `if (iorder=="NNLO")` block:

```fortran
C Find R_Ai
         R_A1=1.0D0
         R_A2=1.0D0
         R_A4=1.0D0
         R_A0=1.0D0
         R_A3=1.0D0
         if (iorder=="NNLO") then
           ... interpolate real R_Ai values from the grid file ...
         endif
```
(`get_yk_new.f:229-305`)

So when `iorder=NLO`, whatever is in the R_Ai file is never read into `R_A0..R_A4` — the defaults of `1.0` (i.e. "no correction") are used regardless. The catch is that the program unconditionally opens and parses the R_Ai file at startup (`get_yk_new.f:116-128`, via `CheckSum`/`PreReadIn`/`ReadIn`) before it even looks at `iorder`, so you can't omit the argument or point it at a missing/empty file — it has to be a well-formed grid file, just not a physically meaningful one.

This isn't a hack invented for this workaround — the original code already has a precedent for it: for W± production, `get_yk_new.f:259-267` hard-wires the interpolation order specifically for the case of "using a dummy R_Ai file in which all the entries are 1."

### R_Ai file format

```
ECM, TYPE_V, PDF
<ECM> <TYPE_V> <PDF_SET>
  Q,qT,y R_A0 R_A1 R_A2 R_A3
<Q> <qT> <y> <R_A0> <R_A1> <R_A2> <R_A3>
<Q> <qT> <y> <R_A0> <R_A1> <R_A2> <R_A3>
...
```

- **Header is exactly 3 lines.** `get_yk_new.f`'s `iHeadLen()` (get_yk_new.f:1369) finds the end of the header dynamically: it reads lines until it finds one whose first 8 characters are `  Q,qT,y` (two leading spaces), and that line counts as part of the header. The two lines before it (`ECM, TYPE_V, PDF` and the actual `<ECM> <TYPE> <PDF>` values) can contain anything parseable — `ReadIn` (get_yk_new.f:530-541) reads the first line and discards it, then does `Read(2,*) ECMC, dummy, dummy` on the second line, so it needs at least one leading numeric field.
- **7 data columns per row**: `Q qT y R_A0 R_A1 R_A2 R_A3`, whitespace-separated. The column count is auto-detected (`Column` subroutine, get_yk_new.f:1458) by scanning until the first column's value (`Q`) repeats — so **consecutive rows must share the same `Q`** (and ideally the same `y`, varying `qT` fastest) the same way the real `R_Ai_*.txt` files and `w_pert_*.out`/`w_asym_*.out` files do. Real Q-major, y-mid, qT-minor grid ordering:
  ```
  91 8.14308 -4.2 1.40419 1.3825 1.18888 1
  91 8.34894 -4.2 1.38738 1.35099 1.17402 1
  ...
  ```
- **Number of rows / grid coverage**: for NLO you don't need the grid to physically match anything — you just need `PreReadIn`/`ReadIn` to successfully parse a rectangular Q/qT/y grid without erroring. The simplest way to guarantee that is to reuse the exact `Q,qT,y` triplets from your `w_pert_*.out` file (same row count, same ordering), so the shapes trivially line up.

### Generating a dummy R_Ai file

`make_dummy_rai.py` (in this repo) builds one automatically from an existing `w_pert_*.out` file's `Q,qT,y` grid, filling in `R_A0..R_A3 = 1` for every row:

```bash
python3 make_dummy_rai.py <w_pert_file> <output_R_Ai_file> [ECM] [TYPE_V] [PDF_SET]

# Example (already generated in this repo as R_Ai_dummy_ZU.txt):
python3 make_dummy_rai.py w_pert_ZU.out R_Ai_dummy_ZU.txt 8000.0 DY CT14nn.00

# Then run NLO as usual:
./get_yk_new w_pert_ZU.out w_asym_ZU.out legacy_y_ZU.out R_Ai_dummy_ZU.txt output_ykR_ZU.out.NLO ZU NLO CT14nn.00
```

For your own custom grid (no matching test data in this repo), point the script at your own `w_pert_*.out` file instead.

### NNLO caveat

If you need to run NNLO without a real R_Ai file, the dummy file above still "works" (the program runs and produces output), but it means **all angular corrections are set to 1, i.e. disabled** — this is a physics choice, not a technical workaround, and will not match a proper NNLO+R_Ai calculation. Only do this if your professor explicitly wants the no-angular-correction NNLO result; otherwise stick to NLO or get a real R_Ai grid from FCFM.

## File Organization

### Essential Files to Copy for Analysis

#### For Drell-Yan (ZU/ZD) at LHC:
```
get_yk_new/
├── get_yk_new            (compiled executable)
├── w_pert_ZU.out         (perturbative input)
├── w_asym_ZU.out         (asymptotic input)
├── legacy_y_ZU.out       (baseline calculation)
└── R_Ai_LHC8_66Q116.txt  (angular functions)
```

#### For W Production at Tevatron:
```
get_yk_new/
├── get_yk_new            (compiled executable)
├── w_pert_wp_tev2.out    (perturbative input)
├── w_asym_wp_tev2.out    (asymptotic input)
├── legacy_y_wp_tev2.out  (baseline calculation)
└── R_Ai_tev2_wp.txt      (angular functions)
```

### Repository Structure

```
get_yk_new/
├── README.md                          # This file
├── CLAUDE.md                          # Development guide
├── Makefile                           # Build configuration
├── make_dummy_rai.py                  # Generates a placeholder R_Ai file for NLO runs
├── get_yk_new.f                       # Main source code
├── get_yk_new_NoRAI                   # Variant without angular corrections
├── get_yk_new_origin                  # Original version
│
├── grids_test_ykR/                    # LHC 8 TeV Drell-Yan test data
│   ├── w_pert_ZU.out
│   ├── w_asym_ZU.out
│   ├── legacy_y_ZU.out
│   ├── R_Ai_LHC8_66Q116.txt
│   ├── legacy_ykR_ZU.out.NLO          # Reference output (NLO)
│   └── legacy_ykR_ZU.out.NNLO         # Reference output (NNLO)
│
├── grids_test_wp_tev2/                # Tevatron W+ production test data
│   ├── w_pert_wp_tev2.out
│   ├── w_asym_wp_tev2.out
│   ├── legacy_y_wp_tev2.out
│   ├── R_Ai_tev2_wp.txt
│   ├── legacy_ykR_wp_tev2.out.NLO
│   └── legacy_ykR_wp_tev2.out.NNLO
│
├── check_R_Ai/                        # Angular function verification
│   └── (verification test case)
│
└── keep_codes/                        # Historical versions
    └── (timestamped backups)
```

## Troubleshooting

### Issue: `make: gfortran: command not found`
**Solution**: Load the GCC module on HPCC
```bash
module load GCC/11.2.0
```

### Issue: `LHAPDF not found`
**Solution**: Verify library path in ~/.bashrc
```bash
echo $LD_LIBRARY_PATH
# Should contain LHAPDF/lib path
lhapdf --version
```

### Issue: `./get_yk_new: cannot open input file`
**Solution**: Verify input files are in the working directory and readable
```bash
ls -lh w_pert*.out w_asym*.out legacy_y*.out R_Ai*.txt
```

### Issue: `Kinematical parameters do not match at line N`
**Cause**: `w_pert`, `w_asym`, and `legacy_y` were generated over different `(Q, qT, y)` grids. See [Grid consistency between w_pert, w_asym, and legacy_y](#grid-consistency-between-w_pert-w_asym-and-legacy_y) for how to diagnose and fix this — usually a mismatched "grid file name" or "Active setting" line in one of the generators' `.in` config files, not the grid definition files (`q_grid.inp`/`qt_grid.inp`/`y_grid.inp`) themselves.

### Issue: Output file contains NaN or -999 values
**Possible Causes**:
- Missing angular correction factors (R_Ai values)
- Mismatched grid points between input files
- Incorrect boson type or order parameter

**Solution**: 
1. Verify grid compatibility: `wc -l w_pert_*.out w_asym_*.out legacy_y_*.out`
2. Check against reference output from `grids_test_ykR/` or `grids_test_wp_tev2/`
3. Review input file headers for parameter matches

### Issue: Performance is slow
**Solution**: Adjust OpenMP threads in ~/.bashrc
```bash
# For better performance on allocated nodes
export OMP_NUM_THREADS=28  # Match --cpus-per-task in SLURM
export OMP_STACKSIZE=2G
```

### Issue: Compilation warnings about unused variables
**Note**: These are safe to ignore in legacy Fortran code. The code has been validated against reference outputs.

## Additional Resources

- **HPCC Documentation**: https://wiki.hpcc.msu.edu/
- **LHAPDF Manual**: https://lhapdf.hepforge.org/
- **MCFM Homepage**: http://mcfm.fnal.gov/
- **ResBos2 Documentation**: Available in `/mnt/home/lopezels/InstallSources/ResBos2/`
- **ROOT Framework**: https://root.cern/

## Support & Contact

For issues specific to:
- **HPCC environment**: Contact MSU HPCC team (hpcc@msu.edu)
- **Code logic/physics**: Review comments in `get_yk_new.f` or consult physics literature on K-factors
- **Input file generation**: See parent research code (MCFM, ResBos, FCFM)

## Citation

If you use this code in research, please cite relevant papers on:
- K-factor calculations and resummation
- The specific physics process (DY, W production)
- PDF sets used in the analysis

---

**Last Updated**: July 2026
**Maintained in**: `/mnt/home/lopezels/Insync/...`
