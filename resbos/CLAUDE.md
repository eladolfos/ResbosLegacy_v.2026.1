# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**ResBos** is a Monte Carlo event generator for vector boson production in hadron collisions (W, Z, photon). It implements CSS resummation formalism with exact matrix elements and spin-correlation effects. The generator produces events at fully differential level, including decay products and quantum correlations.

- **Language**: Fortran 77/90 with C++ interface for ROOT integration
- **Dependencies**: gfortran, g++, ROOT (optional), CERNLIB/LHAPDF (not required in current version)
- **Input**: Grid files (from Legacy code) + configuration file (resbos.in)
- **Output**: ROOT ntuples, HBOOK files, or text histograms

## Architecture

### Core Modules

The code is organized around `resbosm.f` which contains all calculation modules. Key entry points:

- **resbos_root.f** (default): Main program for general vector boson production
- **resbos_root_WZ.f**: Specialized for W/Z production at hadron colliders
- **resbos_root_WZ_cms_zp.f**: W/Z with CMS-specific cuts and Z' handling

### Main Fortran Files

- **resbosm.f** (~658KB): Core MC calculation engine; contains GENERATOR, VEGAS integration, phase space sampling, and common blocks
- **resbosma.f**: Additional subroutines supporting the main generator
- **hvvresph.f**: Handles vector boson (W/Z/photon) resonance production calculations
- **pentagon.f**: Diphoton production amplitudes at loop level (q-qbar-gg → γγ + parton)
- **pncommon.f**: Common blocks for pentagon calculations
- **angular.cc**: Angular distribution calculations (C++)

### Interface Layer

- **froot.c**: C++/ROOT interface providing Fortran-callable functions for ROOT ntuple I/O
  - Functions: `initrootnt_`, `addntbranch_`, `fillntbranch_`, `rootntoutp_`
  - Compiled conditionally with `-DMYROOT` flag (see Makefile lines 40-48)

## Building

```bash
# Build main executable for vector boson production
make resbos_root

# Build W/Z specialized version
make resbos_root_WZ

# Build CMS-specific version with Z' handling
make resbos_root_WZ_cms_zp

# Clean object files and compiled modules
make clean
```

Key compiler flags in Makefile:
- `FC = gfortran` (line 12): Fortran compiler
- `FFLAGS = -O3 -g -fno-automatic` (line 18): Optimization with debugging, static variables disabled
- `CXXFLAGS = -g -std=c++17 -DMYROOT` (line 48): C++ flags with ROOT support

To disable ROOT support and compile without CERN libraries, comment out ROOTDIR, DMYROOT, ROOTLIBS, and ROOTINCLUDE in Makefile (lines 40-43).

## Running

```bash
# Default: reads resbos.in, outputs events
./resbos_root

# Configure via input file: resbos.in
```

### Input Configuration (resbos.in format)

Key parameters in order:
1. **Vegas parameters** (line 1): `nCall, ITMX1, NCALL1, ITMX2, NCALL2, ISEED`
   - ITMX1/NCALL1: preliminary iterations for Vegas sampling
   - ITMX2/NCALL2: main MC iterations (keep ratio ~1000)
   - ISEED: random seed

2. **Grid files** (lines 2-3):
   - Main grid: resummed/NLO portions of d(σ)/d(QT)/dy/dQ²
   - Y-piece grid: high-QT region (optional, use "-" to skip)

3. **Weighting & corrections** (line 4): Un/reweight, K-factors (pT & rapidity), matching flags

4. **Kinematic cuts** (lines 5-7):
   - Line 5: lepton pT, rapidity, ΔR separation
   - Line 6: boson mass range, qT range, boson rapidity
   - Line 7: transverse mass bounds, missing ET

5. **Physics parameters** (line 8): Luminosity, top mass, Higgs mass

6. **Output format** (line 9): ROOTNT1 (ROOT), PAWNT, HBOOK, or GBOOK

7. **Process parameters** (lines 10-14): qT separator, PDF ID, process code, scale choice, boson decay mode, masses/widths, EW flags

## Code Flow

1. **GENERATOR** (in resbosm.f): Sets up integration
2. **Vegas routine**: Samples phase space adaptively
3. **Event generation**:
   - Initialize random number generator with ISEED
   - For each Vegas iteration:
     - Sample 4-momenta of decay products and resummation variables
     - Compute matrix elements (tree-level + NLO corrections)
     - Apply kinematic cuts
     - Fill output histograms/ntuples
4. **Output**: Event distributions to ROOT or HBOOK files

## Key Data Structures (COMMON blocks)

Important common blocks in resbosm.f (lines 71-100):

| Common Block | Purpose |
|---|---|
| `/AA/` | W mass, QCD scales, number of events, orders |
| `/OUT_MOMENTA/` | Final state 4-momenta (lines 76) |
| `/BVEG1/`, `/BVEG2/` | Vegas integration state (lines 77-78) |
| `/STAND0/`, `/STAND1/` | W/Z masses, widths, SM parameters |
| `/MY_SUSY/` | SUSY Higgs parameters (mass, width for H0, H±) |
| `/GridName/`, `/Boson/` | Input grid files and production process |

## Important Implementation Notes

### Warnings in Code

- Hard-wired cuts in WRESPH subroutine optimized for W+ at 1.8 TeV Tevatron (resbos_root.f:27-29)
- Grid files must have consistent PDF sets between main and Y-piece grids (resbos.in:23)
- No automatic correction factor for multiple generator runs (resbos_root.f:122)

### Electroweak Parameters

Four schemes supported (line 26-32 of resbos.in):
- EW_FLAG=0: Original Resbos parameters
- EW_FLAG=1: ZFITTER scheme
- EW_FLAG=2: Effective Born approximation
- EW_FLAG=3,4: FEWZ schemes (complex mass, Gμ scheme)

QED coupling (qed_FLAG): Running α vs. G_μ scheme

CKM matrix: Use FLAG=1 for 2016 PDG values

### Grid File Dependencies

- Grids generated by **Legacy** code (separate utility)
- Required for resummed contributions
- Y-piece grids needed only for high-qT precision (qT > qT_Sep, line 111 of resbos.in)

## Common Development Tasks

**Testing changes to matrix elements**:
- Modify pentagon.f (diphoton) or hvvresph.f (vector bosons)
- Recompile: `make clean && make resbos_root`
- Create small test run with NCALL2=1000, ITMX2=2 to verify fast

**Adding new processes**:
- Create new resbos_root_PROCESSNAME.f variant
- Add compilation target to Makefile
- Ensure grid files available in ./Resbos_grids/

**Debugging numerical issues**:
- Increase verbosity: set IVEGOUT=1 (line 80 in resbosm.f, controlled via resbos.in)
- Check Vegas convergence: NPRN parameter in /BVEG1/ common block
- Validate grid files haven't corrupted (check file size, first line)

## Conventions

- **Masses & widths**: GeV units (resbos.in line 97)
- **Cross sections**: pb (inverse picobarns)
- **Rapidities**: pseudo-rapidity η for leptons/jets
- **4-momenta storage**: (E, px, py, pz) in FINAL_P array, indexed by particle type in NP_TYPE
- **Coupling constants**: α_em at MZ, running α_s from grid files

## External Resources

- Legacy code documentation: http://hep.pa.msu.edu/wwwlegacy/ (grid generation)
- Original ResBos: http://www.pa.msu.edu/~balazs/ResBos
- PAPAGENO base MC structure (resbos_root.f:43-44)
