# w_asym — W-Boson Charge Asymmetry Calculator

Computes the W-boson (and Z-boson) charge asymmetry in hadronic collisions
using CSS/Arnold-Kauffman small-pT resummation expanded to second order in
alpha_s. Supports W+, W-, and Z production for both pp and ppbar colliders.

Original code by Pavel Nadolsky; extended for LHAPDF interface and
kinematic-correction options (see comments in `w_asym.f` and `asym.f`).

---

## Getting the code

```bash
git clone --recursive https://github.com/eladolfos/w_asym_08112022.git
cd w_asym_08112022
```

The `--recursive` flag ensures any submodules are also fetched.

---

## Physics output

For each point on a (Q, qT, y) grid the program writes four columns:

| Column | Quantity |
|--------|----------|
| `asymL0_1` | L_0 piece, particle-1 beam |
| `asymA3_1` | A_3 piece, particle-1 beam |
| `asymL0_1+asymL0_2` | Combined L_0 asymmetry |
| `asymA3_1+asymA3_2` | Combined A_3 asymmetry |

---

## Dependencies

| Library | Version tested | Purpose |
|---------|---------------|---------|
| `gfortran` | system | Fortran compiler |
| `g++` | system | C++ wrapper for LHAPDF |
| LHAPDF 6.x | 6.x | Modern PDF sets (CT18NNLO, etc.) |

The built-in CTEQ packages (`EvlPac02b`, `PrzPac02b`, `QcdPac02b`,
`UtlPac02b`, `EwkPac02b`) and the CT14 Fortran interface (`CT14Pdf.f`) are
already included in the source tree.

---

## Installation on MSU HPCC

### 1. Directory layout

Following the conventions already used on HPCC:

```
/mnt/home/lopezels/
├── InstallSources/
│   └── LHAPDF/          ← LHAPDF installation (headers + lib)
└── SourceCodes/
    └── w_asym/          ← this repository
```

### 2. Install LHAPDF (if not already present)

```bash
cd /mnt/home/lopezels/SourceCodes

# Download LHAPDF 6 from https://lhapdf.hepforge.org
tar xzf LHAPDF-6.X.Y.tar.gz
cd LHAPDF-6.X.Y

./configure --prefix=/mnt/home/lopezels/InstallSources/LHAPDF
make -j4
make install
```

Then make sure these lines are in your `~/.bashrc` (they already appear in
the provided configuration):

```bash
export LD_LIBRARY_PATH=/mnt/home/lopezels/InstallSources/LHAPDF/lib:$LD_LIBRARY_PATH
export PATH=/mnt/home/lopezels/InstallSources/LHAPDF/bin:$PATH
```

Reload your shell:

```bash
source ~/.bashrc
```

Verify:

```bash
lhapdf-config --version
```

### 3. Download a PDF set

```bash
lhapdf install CT18NNLO
```

PDF sets are stored by default under `$prefix/share/LHAPDF/`.

### 4. Clone / copy the source code

```bash
cp -r /path/to/w_asym_08112022 /mnt/home/lopezels/SourceCodes/w_asym
cd /mnt/home/lopezels/SourceCodes/w_asym
```

### 5. Compile

```bash
make clean   # removes *.o, *.co, and the w_asym executable from a previous build
make
```

A successful build produces the executable `w_asym` in the same directory.

If `lhapdf-config` is not on your `PATH`, the Makefile will fail. Confirm
the binary is found:

```bash
which lhapdf-config
```

---

## Running the code

```bash
./w_asym w_asym       # reads w_asym.in, writes w_asym.out
```

Or with a custom job name (the program appends `.in` / `.out`):

```bash
./w_asym myjob        # reads myjob.in, writes myjob.out
```

---

## Input files

### `w_asym.in` — main steering file

```
1960,-1,1                                    > ECM,iBeam(-1=ppbar,1=pp), KinCorr(0/1)
1,1                                          > JWTYPE(jwm=-1,jwp=1,jz=2), JZ_TYPE
0,2,1,1,1                                    > iset, norder(1=LO,2=NLO), muR/Q, muF/Q, iscale
lha_CT18NNLO                                 > PDF file name
./inp/q_grid_test.inp                        > Q grid file
./inp/qt_grid.inp                            > qT grid file
./inp/y_grid.inp                             > rapidity grid file
1   50  5   12 132 30  3 3  1                > Active range (iqmin iqmax iqst  iqtmin iqtmax iqtst  iymin iymax iyst)
```

Key parameters:

| Parameter | Meaning |
|-----------|---------|
| `ECM` | Centre-of-mass energy in GeV (e.g. 1960 for Tevatron) |
| `iBeam` | Beam type: `-1` = ppbar, `+1` = pp |
| `KinCorr` | Kinematic correction: `0` = use Q, `1` = use M_T in x_{1,2} |
| `JWTYPE` | Boson: `-1` = W−, `1` = W+, `2` = Z |
| `iset` | PDF set index (0 when using LHAPDF via `lha_` prefix) |
| `norder` | Perturbative order: `1` = LO, `2` = NLO |
| `muR/Q`, `muF/Q` | Renormalization / factorization scale ratios |

### Grid files (`inp/`)

| File | Content |
|------|---------|
| `q_grid_test.inp` | List of invariant-mass Q values (GeV) |
| `qt_grid.inp` | List of transverse-momentum qT values (GeV) |
| `y_grid.inp` | List of rapidity y values |

Only the grid points selected by the index ranges in `w_asym.in` are
computed.

---

## Output

Results are written to `<jobname>.out`. Each line contains:

```
Q   qT   y   asymL0_1   asymA3_1   asymL0_1+asymL0_2   asymA3_1+asymA3_2
```

Example (Tevatron W+):

```
  80.00    0.1501  2.840    741.14    739.03   -61.18   -60.65
```

---

## Starting a new analysis

Copy these files to your new working directory:

```
w_asym            ← compiled executable
w_asym.in         ← steering file (edit this for your run)
inp/
  q_grid_test.inp ← Q grid  (edit or replace)
  qt_grid.inp     ← qT grid (edit or replace)
  y_grid.inp      ← y grid  (edit or replace)
```

You do **not** need to copy the source or recompile unless you change the
Fortran code. The executable is self-contained once linked against LHAPDF.

Typical workflow:

```bash
# 1. Create a new run directory
mkdir /mnt/home/lopezels/SourceCodes/w_asym/runs/myrun
cd    /mnt/home/lopezels/SourceCodes/w_asym/runs/myrun

# 2. Copy the executable and inputs
cp ../../w_asym .
cp ../../w_asym.in myrun.in
cp -r ../../inp .

# 3. Edit myrun.in (ECM, beam type, PDF, grid index ranges, ...)
#    Edit inp/*.inp to define the desired kinematic grids

# 4. Run
./w_asym myrun

# 5. Output is in myrun.out
```

---

## Source file overview

| File | Description |
|------|-------------|
| `w_asym.f` | Main program: reads input, loops over (Q, qT, y), writes output |
| `asym.f` | Core resummation module (Arnold-Kauffman Y-piece, W/Z/H) |
| `gauss.f` | Gaussian numerical integration |
| `util.f` | String utilities (`trmstr`) |
| `CT14Pdf.f` | Built-in CT14 Fortran PDF interface |
| `lhapdf.cpp` | C++ bridge to LHAPDF 6 |
| `EvlPac02b.for` | CTEQ alpha_s evolution |
| `PrzPac02b.for` | CTEQ PDF parametrization |
| `QcdPac02b.for` | QCD utilities |
| `UtlPac02b.for` | General utilities |
| `EwkPac02b.for` | Electroweak couplings |
| `Makefile` | Build rules |
| `w_asym.in` | Default steering file |
| `inp/` | Default kinematic grids |

---

## HPCC job submission (SLURM example)

```bash
#!/bin/bash
#SBATCH --job-name=w_asym
#SBATCH --output=w_asym_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=4G
#SBATCH --time=04:00:00

source ~/.bashrc
cd /mnt/home/lopezels/SourceCodes/w_asym/runs/myrun

export OMP_NUM_THREADS=28
export OMP_STACKSIZE=2G

./w_asym myrun
```

Submit with:

```bash
sbatch run_w_asym.sb
```
