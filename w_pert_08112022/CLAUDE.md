# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Fortran physics simulation code for calculating differential cross-sections of W/Z boson production at hadron colliders. It implements perturbative QCD calculations with resummation techniques for transverse momentum (qT) distributions.

**Key References:**
- P. Arnold and H. Reno, Nucl. Phys. B319 (1989) 37
- R. Ellis et al., Nucl. Phys. B211 (1983) 106
- K. Kajantie and J. Lindfors, Nucl. Phys. B146 (1978) 465

## Build Commands

### Build
```bash
make clean
make
```
This produces the `w_pert` executable. Compiler: `gfortran` (Makefile uses `COMP=gfortran`). Key flags: `-O2 -g -w -fno-automatic`.

### Requirements
- gfortran (or g77)
- g++ (for C++ components)
- LHAPDF library with development files (lhapdf-config must be available in PATH)

### Clean
```bash
make clean
```
Removes all object files (*.o, *.co), executables (w_pert), and core dumps.

## Running the Code

```bash
./w_pert < w_pert.in
```

Output is written to stdout in columnar format. Redirect to capture results:
```bash
./w_pert < w_pert.in > output.txt
```

### Input Configuration
The input file (`w_pert.in`) contains:
- `ECM`: Center-of-mass energy (GeV) and beam type (iBeam=-1 for ppbar, iBeam=1 for pp)
- `JWTYPE`: Boson type (jwm=-1 for W-, jwp=1 for W+, jz=2 for Z, jph=3 for photon)
- `norder`: Calculation order (1 for LO, 2 for NLO)
- `muR/mT, muF/mT`: Renormalization and factorization scale factors
- `PDF file name`: LHAPDF set (e.g., "lha_CT14nnlo")
- Grid file paths: Q, qT, and y grid definitions in `./inp/`

### Test Run
```bash
./w_pert < test_runs/w_pert_w+.in > output.txt
```
Compare against `test_runs/w_pert_w+.out` for validation.

## Code Architecture

### Main Program
- **w_pert.f**: Main driver. Reads input, manages grids, calls physics modules, formats and prints results.
- Global state shared via COMMON blocks: STU (kinematic variables), STRUCT (parton distributions), JCODES (boson type codes), WPARM/WPARM2 (W/Z mass, couplings).

### Physics Modules
- **pert.f**: Perturbative W/Z production (WXINT module v3.1). Computes dσ/dy/dqT² for given y and qT via parton-level integrals. Returns results in pb/GeV².
  - Top-level routines: `DYPT1` (1st order), `DYPT2` (2nd order), `YMAX` (kinematically allowed rapidity)
  - Uses variable transformation (x,s2) → (zx,zs2 ∈ [0,1]) to smooth integrals
- **gauss.f**: Gaussian numerical integration (2D quadrature)
- **util.f**: Utility functions

### Supporting Modules
- **CT14Pdf.f**: Direct CT14 PDF access
- **lhapdf.cpp**: C++ wrapper to LHAPDF library for PDF member sets
- **Pac files** (.for sources, .o binaries):
  - PrzPac02b: Process packet (hard-scattering kernels)
  - EvlPac02b: Evolution (DGLAP evolution for PDFs)
  - QcdPac02b: QCD (alpha_s, running coupling)
  - UtlPac02b: Utilities
  - EwkPac02b: Electroweak (W/Z masses, couplings, mixing angles)

### Data Files
- `pdf00.pds`, `CT14nnlo.LHgrid`: Parton distribution function data
- `inp/q_grid.inp`, `inp/qt_grid.inp`, `inp/y_grid.inp`: Kinematic grid definitions
- `test_runs/`: Reference output for validation

## Key Variables and COMMON Blocks

### STU block
`ss, tt, uu, qq, qm, qmm, qmu`: Mandelstam variables and invariant masses (GeV²)

### STRUCT block
Parton distributions: `u1, d1, s1, c1, b1, g1` (proton) and `u2, d2, ...` (antiproton/proton)

### WPARM, WPARM2 blocks
W/Z parameters: masses (WMASS, ZMASS), couplings (GGU, GGD, GVU, GVD, GAU, GAD), weak mixing angle sin²θ_W

### JCODES block
Boson type flags: jwm=-1 (W⁻), jwp=1 (W⁺), jz=2 (Z), jph=3 (γ)

## Precision and Accuracy Notes

- **Integration tolerance**: `ERREL` (default 0.002 to 0.005) controls relative error in Gaussian quadrature. Set in w_pert.f COMMON /ERR/
- **Scale choices**: Renormalization/factorization scales default to qT mass (iscale=1) or Q² (iscale=0)
- **Resummation matching**: Perturbative and resummed calculations must cancel to high precision at intermediate qT; this requires consistent structure function approximations

## Compiler Flags in Makefile

- `-O2`: Optimization level
- `-g`: Include debugging symbols
- `-w`: Suppress compiler warnings (legacy code)
- `-fno-automatic`: Static allocation of local variables (Fortran 77 style)

Adjust these flags in `COMPFLAGS` for debugging, profiling, or alternative optimization levels.

## Common Modifications

1. **Change PDF set**: Edit `w_pert.in`, line 4 (e.g., "lha_CT14lo" for LO, "lha_CT14nnlo" for NNLO)
2. **Adjust beam type**: Edit `w_pert.in`, line 1 (iBeam=-1 for ppbar, iBeam=1 for pp)
3. **Modify scales**: Edit `w_pert.in`, line 3 (muR/mT, muF/mT factors, iscale)
4. **Change integration tolerance**: Edit w_pert.f, COMMON /ERR/ `ERREL` value
5. **Output format**: Controlled by format statements in main loop (w_pert.f). Currently: Q, qT, y, pert_1, pert_2, K-factor

## Debugging

- **Recompile with debug symbols**: Set `COMPFLAGS=-c -g` in Makefile (current default includes `-g`)
- **Suppress warnings**: Makefile uses `-w` flag; remove to see warnings
- **Print intermediate values**: Add write statements to pert.f or w_pert.f near calculations of interest
- **Numerical instabilities**: Check ERREL tolerance; reduce for problematic kinematics
- **Segmentation faults**: May indicate array bounds violations; check grid dimensions (inp/ files) vs. array sizes in COMMON blocks
