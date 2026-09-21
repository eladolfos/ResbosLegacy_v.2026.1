# w_pert — W/Z Boson Differential Cross-Section Calculator

Computes the differential cross-sections (dσ/dy/dqT²) for W and Z boson
production in hadronic collisions using CSS (Collinear-Soft-Soft) resummation
expanded to second order in alpha_s. Supports W±, Z production for both pp
and ppbar colliders (Tevatron and LHC).

Original code by Pavel Nadolsky; extended for LHAPDF interface and kinematic
corrections (see comments in `w_pert.f` and `pert.f`).

---

## Physics output

For each point on a (Q, qT, y) grid the program writes columns:

| Column | Quantity |
|--------|----------|
| `Q` | Invariant mass of boson (GeV) |
| `qT` | Transverse momentum (GeV) |
| `y` | Rapidity |
| `pert_1` | 1st-order (LO) contribution |
| `pert_2` | 2nd-order (NLO) contribution |
| `factor_K` | K-factor = (pert_1 + pert_2) / pert_1 |

Results in pb/GeV².

---

## Dependencies

| Library | Version tested | Purpose |
|---------|---------------|---------|
| `gfortran` | system (HPCC default) | Fortran 77/90 compiler |
| `g++` | system (HPCC default) | C++ wrapper for LHAPDF |
| LHAPDF 6.x | 6.x+ | Modern PDF sets (CT14NNLO, CT18NNLO, etc.) |

The built-in CTEQ packages (`EvlPac02b`, `PrzPac02b`, `QcdPac02b`,
`UtlPac02b`, `EwkPac02b`) and the CT14 Fortran interface (`CT14Pdf.f`) are
already included in the source tree.

---

## Installation on MSU HPCC

### 1. Directory layout

Following the conventions already established on your HPCC account:

```
/mnt/home/lopezels/
├── InstallSources/
│   ├── LHAPDF/              ← LHAPDF installation (headers + lib)
│   ├── HOPPET1/             ← (existing)
│   ├── ROOT/                ← (existing)
│   ├── ResBos2/             ← (existing)
│   └── MCFM-10.3/           ← (existing)
└── SourceCodes/
    ├── w_asym_08112022/     ← (existing)
    └── w_pert_08112022/     ← this repository (copy here)
```

### 2. Verify LHAPDF installation

Your `.bashrc` already contains LHAPDF paths. Verify it is installed:

```bash
source ~/.bashrc
which lhapdf-config
lhapdf-config --version
```

If LHAPDF is not installed, follow the installation steps in the w_asym README.

### 3. Download/copy the source code

If copying from your local machine:

```bash
cp -r /path/to/w_pert_08112022 /mnt/home/lopezels/SourceCodes/w_pert
cd /mnt/home/lopezels/SourceCodes/w_pert
```

Or if cloning from GitHub (once available):

```bash
cd /mnt/home/lopezels/SourceCodes
git clone https://github.com/yourusername/w_pert_08112022.git
cd w_pert_08112022
```

### 4. Update Makefile for LHAPDF (if needed)

The Makefile queries `lhapdf-config` to find LHAPDF headers and libraries.
On HPCC, after sourcing `.bashrc`, this should work automatically.

To verify the Makefile will find LHAPDF:

```bash
lhapdf-config --incdir
lhapdf-config --ldflags
lhapdf-config --cppflags
```

If any command fails, re-source your `.bashrc`:

```bash
source ~/.bashrc
```

**Optional fix**: If `lhapdf-config` is still not found, edit the Makefile to
explicitly point to LHAPDF:

```makefile
# Replace this line:
LHAPDFINCLUDE = $(shell lhapdf-config --incdir)

# With an absolute path:
LHAPDFINCLUDE = /mnt/home/lopezels/InstallSources/LHAPDF/include
```

And similarly for `LHAPDFLIBS`.

### 5. Compile

```bash
source ~/.bashrc  # Ensure LHAPDF paths are set

cd /mnt/home/lopezels/SourceCodes/w_pert

make clean        # removes old *.o, *.co, and w_pert executable
make              # compile all sources and link
```

A successful build produces the executable `w_pert` in the same directory.

**Troubleshooting**:

- **"lhapdf-config not found"**: Ensure `~/.bashrc` is sourced and LHAPDF is on PATH.
- **Linker errors with LHAPDF**: Add `-lLHAPDF` explicitly in the Makefile `LHAPDFLIBS` line.
- **Compiler not found**: HPCC default modules should provide `gfortran` and `g++`. If not:
  ```bash
  module load icc  # or appropriate compiler module
  ```

---

## Running the code

### Basic execution

```bash
./w_pert < w_pert.in
```

Output is written to stdout. Redirect to save results:

```bash
./w_pert < w_pert.in > w_pert.out
```

### Input file format

The main input file `w_pert.in` controls the calculation:

```
1960,-1                                   > ECM,iBeam(-1=ppbar,1=pp)
1,0,1                                     > JWTYPE(jwm=-1,jwp=1,jz=2,jph=3,JHB=4), IDO_CBAR(0,1), JZ_TYPE(1,-1,0)
0,2,1,1,1                                 > iset, norder(1=LO,2=NLO), muR/mT, muF/mT, iscale(0=Q2,1=MT)
lha_CT14nnlo                              > PDF file name
./inp/q_grid_test.inp                     > Q grid file
./inp/qt_grid.inp                         > qT grid file
./inp/y_grid.inp                          > rapidity grid file
1  105  1   12 132 12  1 7  2             > Active range (iqmin iqmax iqst iqtmin iqtmax iqtst iymin iymax iyst)
```

**Key parameters:**

| Parameter | Meaning |
|-----------|---------|
| `ECM` | Centre-of-mass energy in GeV (e.g., 1960 for Tevatron, 7000 for LHC 7 TeV) |
| `iBeam` | Beam type: `-1` = ppbar (Tevatron), `+1` = pp (LHC) |
| `JWTYPE` | Boson: `-1` = W−, `1` = W+, `2` = Z, `3` = γ, `4` = Higgs |
| `IDO_CBAR` | Heavy-flavor option: `0` = standard, `1` = include c-bar terms |
| `JZ_TYPE` | Z asymmetry: `1` = up-quark, `-1` = down-quark, `0` = combined |
| `iset` | PDF set index (0 when using LHAPDF via `lha_` prefix) |
| `norder` | Perturbative order: `1` = LO, `2` = NLO |
| `muR/mT, muF/mT` | Renormalization / factorization scale ratios (1.0, 0.5, 0.25, 2.0 are common) |
| `iscale` | Scale type: `0` = Q², `1` = M_T (transverse mass) |

### Grid files (`inp/`)

Three kinematic grids define the calculation points:

| File | Content | Example |
|------|---------|---------|
| `q_grid.inp` | List of invariant-mass Q values (GeV) | 60, 80, 90, 100, ... |
| `qt_grid.inp` | List of transverse-momentum qT values (GeV) | 0.1, 0.2, 0.5, 1.0, ... |
| `y_grid.inp` | List of rapidity y values | -3.0, -2.5, ..., 2.5, 3.0 |

Only grid points selected by index ranges in `w_pert.in` (lines 7) are computed.

**Grid file format** (plain text, one value per line):

```
60.0
70.0
80.0
...
```

---

## Output format

Results are written to stdout (or redirected to `.out` file). Each line contains:

```
Q       qT      y       pert_1          pert_2          K_factor        (additional columns)
```

**Example (Tevatron W+ at Q=80, qT=0.5, y=1.5)**:

```
 80.00    0.5000   1.500        125.34      -10.52    0.9161     ...
```

Units: cross-sections in pb/GeV², energies in GeV.

---

## Starting a new analysis

### Setup a new run directory

```bash
# Create directory for your analysis
mkdir -p /mnt/home/lopezels/SourceCodes/w_pert/runs/myrun_20250702
cd /mnt/home/lopezels/SourceCodes/w_pert/runs/myrun_20250702

# Copy the executable and input files
cp ../../w_pert .
cp ../../w_pert.in myrun.in
cp -r ../../inp .
cp ../../test_runs/*.out .  # (optional: reference outputs)
```

### Customize your input

Edit `myrun.in` to set:

- ECM and beam type (iBeam)
- Boson type (JWTYPE)
- PDF set (`lha_CT14nnlo`, `lha_CT18NNLO`, etc.)
- Calculation order (norder = 1 for LO, 2 for NLO)
- Scale factors (muR/mT, muF/mT)
- Grid ranges (if not computing all points)

### Run the calculation

```bash
./w_pert < myrun.in > myrun.out
cat myrun.out | head -20  # Peek at results
```

### Postprocess results

Results are in plain text. Pipe to `awk` or `grep` to extract specific Q/qT/y:

```bash
# Extract results for Q=80 GeV
grep "80.00" myrun.out | head -10

# Extract results for all qT < 1.0 GeV
awk '$2 < 1.0 {print}' myrun.out
```

---

## Files to copy for analysis

**Minimal set for running pre-compiled code:**

```
w_pert              ← compiled executable
w_pert.in           ← steering file (edit for your run)
inp/
  q_grid.inp        ← Q grid
  qt_grid.inp       ← qT grid
  y_grid.inp        ← y grid
```

You do **not** need source files (*.f, *.for, *.cpp) or the Makefile unless
you modify the Fortran code.

**Complete set for development:**

```
w_pert              ← executable (or recompile from source)
w_pert.f            ← main program
pert.f              ← perturbative W/Z module
gauss.f             ← integration
util.f              ← utilities
CT14Pdf.f           ← PDF interface
lhapdf.cpp          ← LHAPDF C++ wrapper
*.for               ← CTEQ packet files (EvlPac02b, PrzPac02b, etc.)
Makefile            ← build rules
w_pert.in           ← steering file
inp/                ← grid files
test_runs/          ← reference outputs
```

---

## HPCC job submission (SLURM example)

Create a file `submit_w_pert.sb`:

```bash
#!/bin/bash
#SBATCH --job-name=w_pert
#SBATCH --output=w_pert_%j.log
#SBATCH --error=w_pert_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1        # w_pert is single-threaded
#SBATCH --mem=2G
#SBATCH --time=04:00:00          # Adjust as needed
#SBATCH --partition=general      # or your preferred partition

# Load environment
source ~/.bashrc

# Ensure LHAPDF is on path
export LD_LIBRARY_PATH=/mnt/home/lopezels/InstallSources/LHAPDF/lib:$LD_LIBRARY_PATH
export PATH=/mnt/home/lopezels/InstallSources/LHAPDF/bin:$PATH

# Run directory
cd /mnt/home/lopezels/SourceCodes/w_pert/runs/myrun

# Execute
./w_pert < myrun.in > myrun.out

echo "w_pert job completed"
```

Submit with:

```bash
sbatch submit_w_pert.sb
```

Monitor:

```bash
squeue -u lopezels
```

---

## Source file overview

| File | Description |
|------|-------------|
| `w_pert.f` | Main program: reads input, manages grids, calls pert.f, formats output |
| `pert.f` | WXINT module v3.1: core resummation (Arnold-Kauffman), W/Z cross-sections |
| `gauss.f` | Gaussian 2D numerical integration (quadrature) |
| `util.f` | Utility functions (string trimming, etc.) |
| `CT14Pdf.f` | Built-in CT14 Fortran PDF interface |
| `lhapdf.cpp` | C++ bridge to LHAPDF 6 |
| `EvlPac02b.for` | CTEQ PDF evolution (DGLAP) |
| `PrzPac02b.for` | CTEQ PDF parametrization |
| `QcdPac02b.for` | QCD utilities (alpha_s) |
| `UtlPac02b.for` | General utility functions |
| `EwkPac02b.for` | Electroweak couplings (W/Z masses, mixing angles) |
| `Makefile` | Build rules |
| `CLAUDE.md` | Developer reference for code architecture |
| `w_pert.in` | Default steering file |
| `inp/` | Default kinematic grids |
| `test_runs/` | Reference input/output for validation |

---

## Validation

Test the installation with the provided reference run:

```bash
cd /mnt/home/lopezels/SourceCodes/w_pert
./w_pert < test_runs/w_pert_w+.in > test_output.out

# Compare first 20 lines against reference
head -20 test_output.out
head -20 test_runs/w_pert_w+.out

# Check for differences (should be negligible, < 1%)
diff test_output.out test_runs/w_pert_w+.out | head
```

---

## Troubleshooting

### Compilation fails with "lhapdf-config not found"

```bash
source ~/.bashrc
which lhapdf-config
make clean && make
```

### Runtime error: "Cannot find PDF file"

Ensure the PDF set is installed:

```bash
lhapdf list | grep CT14nnlo
lhapdf install CT14nnlo   # if not present
```

Update `w_pert.in` line 4 to match an installed set:

```
lha_CT14nnlo         # available
```

### Segmentation fault or NaN results

- Check grid file paths in `w_pert.in` (lines 5-7)
- Verify grid file format: one value per line, no header
- Reduce integration tolerance: edit `w_pert.f`, COMMON /ERR/ `ERREL`
- Check kinematic constraints (e.g., qT should not exceed Q)

### Output format is truncated or wrong

Edit the format statements in `w_pert.f` main loop. Search for `write(*,'(` and
adjust field widths (e.g., change `f8.3` to `f10.4` for higher precision).

---

## References

See `CLAUDE.md` for detailed physics module descriptions, compiler flags, and
key COMMON blocks.

- P. Arnold and H. Reno, Nucl. Phys. B319 (1989) 37
- R. Ellis et al., Nucl. Phys. B211 (1983) 106
- K. Kajantie and J. Lindfors, Nucl. Phys. B146 (1978) 465
