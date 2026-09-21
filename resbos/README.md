# ResBos Installation and Setup Guide for HPCC (Michigan State University)

This guide provides step-by-step instructions for compiling and running ResBos on Michigan State University's High Performance Computing Cluster (HPCC).

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [HPCC Environment Setup](#hpcc-environment-setup)
3. [Installation Steps](#installation-steps)
4. [Configuration](#configuration)
5. [Running on HPCC](#running-on-hpcc)
6. [Troubleshooting](#troubleshooting)

---

## Prerequisites

ResBos requires the following tools to be installed and available on HPCC:

### Required
- **gfortran** (Fortran 77/90 compiler) — provided by module `GCC`
- **g++** (C++ compiler) — provided by module `GCC`
- **GNU Make** — typically available by default

### Optional but Recommended
- **ROOT** (for ntuple output) — available as module or via custom installation
- **LHAPDF** (PDF library) — custom installation recommended for advanced features

---

## HPCC Environment Setup

### Step 1: Load Required Modules

Add the following to your `~/.bashrc` file on HPCC:

```bash
# Load compiler module
module load GCC/11.2.0

# Load ROOT (choose appropriate version)
# Option A: If ROOT module is available on HPCC
module load ROOT/6.26.06-GCCcore-11.2.0

# Option B: If using custom ROOT installation
source /mnt/home/YOUR_USERNAME/InstallSources/ROOT/bin/thisroot.sh
```

### Step 2: Add Optional Physics Libraries (if available)

If you have custom installations of LHAPDF, HOPPET, or MCFM, add to `~/.bashrc`:

```bash
# LHAPDF configuration
export LD_LIBRARY_PATH=/mnt/home/YOUR_USERNAME/InstallSources/LHAPDF/lib:$LD_LIBRARY_PATH
export PATH=/mnt/home/YOUR_USERNAME/InstallSources/LHAPDF/bin:$PATH

# HOPPET configuration (if using for PDF evolution)
export LD_LIBRARY_PATH=/mnt/home/YOUR_USERNAME/InstallSources/HOPPET1/lib:$LD_LIBRARY_PATH
export PATH=/mnt/home/YOUR_USERNAME/InstallSources/HOPPET1/bin:$PATH
```

### Step 3: Load Bashrc and Verify Setup

```bash
source ~/.bashrc
which gfortran      # Should show path to GCC gfortran
which g++           # Should show path to GCC g++
root-config --version   # Should show ROOT version (if installed)
```

---

## Installation Steps

### Step 1: Create Installation Directory

```bash
# Create a directory for ResBos (adjust path as needed)
mkdir -p /mnt/home/YOUR_USERNAME/InstallSources/ResBos
cd /mnt/home/YOUR_USERNAME/InstallSources/ResBos
```

### Step 2: Copy ResBos Source Code

Copy the ResBos source files to your HPCC home directory:

```bash
# If not already there, copy all .f, .c, and Makefile
cp /path/to/resbos/*.f .
cp /path/to/resbos/*.c .
cp /path/to/resbos/Makefile .
cp /path/to/resbos/*.in .
```

### Step 3: Create Grid Directory

```bash
mkdir -p Resbos_grids
# Download or copy grid files from Legacy code output
# Grid files should be placed in ./Resbos_grids/
```

### Step 4: Configure Makefile for HPCC

Edit the `Makefile` to ensure proper compiler and library settings:

```makefile
# Verify these lines (around line 12-18):
FC = gfortran          # Fortran compiler
CXX = g++              # C++ compiler
CXXFLAGS0 = -g -O3     # Compiler optimization level
FFLAGS = -O3 -g -fno-automatic

# For ROOT support, verify (around line 40-43):
ROOTDIR = $(shell root-config --prefix)
DMYROOT = -DMYROOT
ROOTLIBS = $(shell root-config --libs) -lstdc++
ROOTINCLUDE = -I $(shell root-config --incdir)
```

**Note:** If ROOT is not available, comment out the ROOTDIR, DMYROOT, ROOTLIBS, and ROOTINCLUDE lines.

### Step 5: Compile ResBos

```bash
# Clean any previous builds
make clean

# Build main executable
make resbos_root

# Alternative: build specialized versions
# For W/Z production:
make resbos_root_WZ

# For CMS analysis with Z':
make resbos_root_WZ_cms_zp

# Verify executable exists
ls -lh resbos_root
```

### Step 6: Test Compilation

Run a quick test to verify the executable works:

```bash
# Create a test input file (resbos.in in current directory)
./resbos_root < /dev/null
# Should not crash with compilation errors
```

---

## Configuration

### Input Configuration File: `resbos.in`

ResBos reads input parameters from `resbos.in`. Key parameters to configure:

```fortran
! Line 1: Vegas integration parameters
1,30,30000,50,500000,182019    ! nCall, ITMX1, NCALL1, ITMX2, NCALL2, ISEED

! Lines 2-3: Grid files (adjust path to your Resbos_grids directory)
./Resbos_grids/legacy_zd_w321_8TeV.out       ! Main grid
./Resbos_grids/legacy_zd_ykR_8TeV.out        ! Y-piece grid

! Line 4: Weighting and corrections
0, 0, 0, 2, 1, 0              ! Reweight, K-factors, matching

! Lines 5-7: Kinematic cuts
0.0, -10, 0.0, 10000.0, 10    ! pT, rapidity, separation cuts
66.0, 116.0, 0.0, 50., -99.0, 99.0  ! Mass and qT cuts
0.0, 10000.,0.0               ! Transverse mass cuts

! Line 8: Physics parameters
1.0, 173.5, 125.0             ! Luminosity (pb^-1), mt, mH

! Line 9: Output format
ROOTNT1                        ! ROOT output (requires ROOT)

! Line 10: Process parameters
1.5, 903, 0, 5, 0, 0, 0, 0.0  ! qT_Sep, PDF ID, Process, Scale, etc.

! Line 11: Boson decay
WW,0                           ! Decay mode, additional flag

! Lines 12-13: Boson masses and widths (PDG 2023 values)
80.385, 2.0906                 ! W mass (GeV), width
91.1876, 2.4952                ! Z mass (GeV), width

! Line 14: Electroweak scheme
1, 1, 0, 0                     ! EW_FLAG, qed_FLAG, CKM_FLAG, KFAC_ANGFUNC
```

**Important:** Grid files must exist in the `./Resbos_grids/` directory. These are generated by the Legacy code.

### Adjusting for Your Analysis

- **Luminosity** (line 8, first value): Set to your experiment's integrated luminosity
- **Process** (line 10, third value): See resbos_root.f documentation for process codes
- **PDF ID** (line 10, second value): 903 = CTEQ6, or use LHAPDF ID if configured
- **Cuts** (lines 5-7): Adjust to match your analysis requirements

---

## Running on HPCC

### Interactive Run (Small Tests)

```bash
# Submit interactive job
qsub -I -l nodes=1:ppn=28 -l walltime=01:00:00 -q short

# Load environment and run
module load GCC/11.2.0
source /path/to/thisroot.sh
cd /mnt/home/YOUR_USERNAME/InstallSources/ResBos

# Run with default resbos.in
./resbos_root

# Run with custom input file
./resbos_root < custom_input.in
```

### Batch Job Submission

Create a PBS/Torque job script (`submit_resbos.pbs`):

```bash
#!/bin/bash
#PBS -N resbos_production
#PBS -l nodes=1:ppn=28
#PBS -l walltime=04:00:00
#PBS -l mem=120gb
#PBS -q main
#PBS -j oe
#PBS -o resbos_${PBS_JOBID}.log

# Load modules
module load GCC/11.2.0
source /path/to/thisroot.sh

# Set OpenMP threads
export OMP_NUM_THREADS=28
export OMP_STACKSIZE=2G

# Navigate to work directory
cd /mnt/home/YOUR_USERNAME/InstallSources/ResBos

# Run ResBos
./resbos_root > output_${PBS_JOBID}.txt 2>&1

# Optional: Copy output to storage
# cp *.root /mnt/gs21/scratch/YOUR_USERNAME/
```

### Submit Job

```bash
qsub submit_resbos.pbs

# Check job status
qstat -u YOUR_USERNAME

# Monitor job output (while running)
tail -f resbos_*.log
```

### Parallel Runs (Multiple Grid Points)

For parameter scans, submit multiple jobs with different inputs:

```bash
#!/bin/bash
# create_jobs.sh

for seed in {1..10}; do
    # Create input file
    sed "s/SEED_PLACEHOLDER/$seed/" resbos_template.in > resbos_${seed}.in
    
    # Create job script
    cat > job_${seed}.pbs << EOF
#!/bin/bash
#PBS -N resbos_${seed}
#PBS -l nodes=1:ppn=28
#PBS -l walltime=04:00:00
#PBS -q main
module load GCC/11.2.0
source /path/to/thisroot.sh
cd /mnt/home/YOUR_USERNAME/InstallSources/ResBos
./resbos_root < resbos_${seed}.in > output_${seed}.log 2>&1
EOF
    
    # Submit job
    qsub job_${seed}.pbs
done
```

---

## Output Files

After running, ResBos generates:

- **`*.root`** — ROOT ntuples (if ROOTNT1 output format selected)
- **`*.dat`** — Histogram data files
- **Console output** — Integration convergence info and cross-section summary

### Processing ROOT Output

```bash
# Load ROOT and analyze output
root -l output_file.root

# In ROOT shell:
root [0] TBrowser b;  // Opens ROOT file browser
root [0] events->Scan();  // View ntuple contents
```

---

## Troubleshooting

### Compilation Errors

#### "gfortran: command not found"
```bash
# Load GCC module
module load GCC/11.2.0
module list  # Verify GCC is loaded
```

#### "root-config: command not found"
```bash
# Option 1: Load ROOT module
module load ROOT/6.26.06-GCCcore-11.2.0

# Option 2: Use custom ROOT installation
source /mnt/home/YOUR_USERNAME/InstallSources/ROOT/bin/thisroot.sh
```

#### "undefined reference to `TFile'" or ROOT linking errors
```bash
# Ensure ROOT was compiled with same GCC version
# Verify in Makefile:
$(shell root-config --libs)      # Should include -lCore, -lTree, etc.
$(shell root-config --cflags)    # Should show correct include paths

# Recompile with verbose output
make resbos_root 2>&1 | tail -50
```

### Runtime Errors

#### "Grid file not found"
```bash
# Verify grid files exist and are readable
ls -lh Resbos_grids/
file Resbos_grids/*.out  # Check file type

# Update paths in resbos.in to absolute or relative paths
```

#### "Permission denied" on HPCC
```bash
# Ensure files have correct permissions
chmod +x resbos_root
chmod +r resbos.in
chmod -R +r Resbos_grids/
```

#### "Segmentation fault" during execution
```bash
# Increase stack size on HPCC
export OMP_STACKSIZE=4G  # Increase from 2G

# Recompile with debugging symbols
make clean
make resbos_root  # Makefile uses -g flag by default

# Run with gdb (if available)
gdb ./resbos_root
(gdb) run < resbos.in
(gdb) backtrace  # Show call stack at crash
```

### Performance Issues

#### Slow Vegas convergence
```bash
# Reduce precision requirements in resbos.in
# Decrease ITMX1, NCALL1, ITMX2 (fewer iterations)
# This trades accuracy for speed (useful for testing)
```

#### Out of memory errors
```bash
# Request more memory in PBS script
#PBS -l mem=240gb  # Increase memory

# Or reduce number of events (NCALL parameters)
```

---

## Environment File Template

For convenience, save the following as `load_resbos_env.sh`:

```bash
#!/bin/bash
# load_resbos_env.sh - Load ResBos environment on HPCC

MODULE_GCC="GCC/11.2.0"
MODULE_ROOT="ROOT/6.26.06-GCCcore-11.2.0"
RESBOS_DIR="/mnt/home/$(whoami)/InstallSources/ResBos"

# Load modules
module load ${MODULE_GCC}
module load ${MODULE_ROOT}

# Alternative: use custom ROOT
# source /mnt/home/$(whoami)/InstallSources/ROOT/bin/thisroot.sh

# Set OpenMP settings for 28-core nodes
export OMP_NUM_THREADS=28
export OMP_STACKSIZE=2G

# Add to PATH
export PATH=${RESBOS_DIR}:${PATH}

echo "ResBos environment loaded successfully"
echo "RESBOS_DIR: ${RESBOS_DIR}"
module list | grep -E "GCC|ROOT"
```

Usage:
```bash
source load_resbos_env.sh
cd ${RESBOS_DIR}
./resbos_root
```

---

## References

- **ResBos Repository**: http://www.pa.msu.edu/~balazs/ResBos
- **Legacy Grid Generation**: http://hep.pa.msu.edu/wwwlegacy/
- **HPCC Documentation**: https://wiki.hpcc.msu.edu/
- **ROOT Documentation**: https://root.cern/
- **CTEQ Collaboration**: http://cteq.org/

---

## Support and Contact

For issues specific to:
- **ResBos code**: See resbos_root.f comments or contact Pavel Nadolsky (nadolsky@pa.msu.edu)
- **HPCC resources**: Contact MSU HPCC support (hpcc-support@msu.edu)
- **Your modifications**: Refer to CLAUDE.md for code architecture
