# couplings.f90 (from src/dish)

## Overview
The `couplings` module is central to handling Nonadiabatic Couplings (NACs) and associated eigenenergies within the 'dish' simulation package. It defines a data structure `overlap` to store these quantities over time. The module provides functionalities to:
1.  Calculate NACs from pairs of VASP WAVECAR files using a finite difference formula.
2.  Read and write coupling data from/to a binary file named `COUPCAR`.
3.  Read and write coupling data from/to formatted text files (`EIGTXT` for energies, `NATXT` or `NATXT_RE`/`NATXT_IM` for couplings).
4.  Initialize and manage WAVECAR file information for NAC calculations.

## Key Components
- **`module couplings`**:
    - Uses modules `prec` (precision), `constants` (physical constants like `imgUnit`), `lattice` (lattice information), `wavecar` (WAVECAR reading utilities like `sysinfo`, `LOADWAVE`, `freemem`, `closewav`), and `fileio` (possibly for `namdInfo` type or other I/O).
- **`type overlap`**:
    - `NBANDS`: Integer, number of bands.
    - `TSTEPS`: Integer, number of time steps.
    - `dt`: Real(q), time step duration.
    - `Dij`: Complex(q), allocatable 3D array (`NBANDS, NBANDS, TSTEPS`), stores coupling matrix D_ij(t).
    - `DijR`, `DijI`: Real(q), allocatable 3D arrays, for real and imaginary parts of `Dij` when read from separate files.
    - `Eig`: Real(q), allocatable 2D array (`NBANDS, TSTEPS`), stores eigenenergies E_i(t).

### Subroutines for Binary COUPCAR I/O
- **`CoupFromFile(olap)`**:
    - Reads coupling data (`Dij`, `Eig`) into an `overlap` type variable `olap` from a binary file named `COUPCAR`.
    - Checks for consistency between `olap` dimensions and metadata stored in `COUPCAR`.
- **`CoupToFile(olap)`**:
    - Writes coupling data from an `overlap` type variable `olap` to `COUPCAR`.
    - Stores metadata (record length, nbands, nsw, dt) in the first record.

### Subroutines for Text File I/O (Energies and NACs)
- **`writeNaEig(olap, inp)`**:
    - Writes energies from `olap%Eig` to `EIGTXT` and complex NACs from `olap%Dij` to `NATXT`.
    - Only writes data for bands specified by `inp%BMIN` to `inp%BMAX`.
- **`readNaEig(olap_sec, inp)`**:
    - Reads energies into `olap_sec%Eig` from `EIGTXT` and complex NACs into `olap_sec%Dij` from `NATXT`.
- **`readNaEigCpx(olap_sec, inp)`**:
    - Reads energies into `olap_sec%Eig` from `EIGTXT`.
    - Reads real part of NACs into `olap_sec%DijR` from `NATXT_RE`.
    - Reads imaginary part of NACs into `olap_sec%DijI` from `NATXT_IM`.
    - Combines them into `olap_sec%Dij = olap_sec%DijR + olap_sec%DijI * imgUnit`.
    - Deallocates `DijR` and `DijI`.

### Subroutines for NAC Calculation from WAVECARs
- **`initAB(FileA, FileB, waveA, waveB)`**:
    - Initializes `waveinfo` types `waveA` and `waveB` for two WAVECAR files (`FileA`, `FileB`).
    - Assigns file units and calls `sysinfo` to read metadata.
    - Performs consistency checks (same ispin, nkpts, nbands, nplws).
- **`CoupIJ(waveA, waveB, Cij)`**:
    - Calculates NACs `Cij` between states from `waveA` (time t) and `waveB` (time t+dt).
    - Formula: `Cij(i,j) = (<psi_i(t)|psi_j(t+dt)> - <psi_j(t)|psi_i(t+dt)>)`. Note: Division by `2*dt` is not done here.
    - Loads all band coefficients from `waveA` into `crA` and from `waveB` into `crB` using `LOADWAVE`.
    - `Cij(j,i) = -DCONJG(Cij(i,j))`.
- **`finishAB(waveA, waveB)`**:
    - Frees memory associated with `waveA`, `waveB` (via `freemem`) and closes their file units (via `closewav`).

### Main Orchestration Subroutine
- **`TDCoupIJ(rundir, olap, olap_sec, inp)`**:
    - Time-Dependent Couplings: Initializes `olap` and `olap_sec` (likely a sub-section of bands from `olap`).
    - Checks if `COUPCAR` exists.
        - If yes, and `inp%LCPTXT` is true, it calls `readNaEig` (or `readNaEigCpx` if `inp%LCPTROT` is true, though `readNaEigCpx` is called from `main.f90` based on typical flow, not directly here based on `LCPTROT`). If `inp%LSPACE` is true, calls `initspace`.
        - If `COUPCAR` exists but `inp%LCPTXT` is false, it indicates an unsupported path (stops, recommending text files).
    - If `COUPCAR` does not exist, it indicates an unsupported path (stops, recommending text files). The code block for calculating from WAVECARs is commented out, suggesting a workflow shift towards pre-calculated text files.
    - Deallocates `olap%Dij` and `olap%Eig` at the end (as `olap_sec` now holds the relevant data).

- **`initspace(olap_sec, inp)`**:
    - Seems to process energies (`olap_sec%Eig`) based on `inp%ACBASIS` (active basis?) and `inp%NELECTRON` to compute new energies for a different basis set (`olap_acbas%Eig`).
    - Writes these new energies to `ACEIGTXT`.

## Important Variables/Constants
- **`COUPCAR`**: Default binary file name for storing/retrieving NACs and energies.
- **`EIGTXT`**: Text file for eigenenergies.
- **`NATXT`**: Text file for complex NACs.
- **`NATXT_RE`, `NATXT_IM`**: Text files for real and imaginary parts of NACs, respectively.
- **`ACEIGTXT`**: Text file for energies in an "active basis" calculated by `initspace`.
- **`waveinfo` type**: Defined in `wavecar` module, holds WAVECAR metadata and planewave coefficients.
- **`namdInfo` type**: Likely defined in `fileio` or another module, holds input parameters like `NBANDS`, `NSW` (number of steps), `POTIM` (time step), `BMIN`, `BMAX` (band range), `LCPTXT` (logical flag for using text coupling files).

## Usage Examples
This module is used internally by the simulation program. `TDCoupIJ` is the main entry point called to load or prepare the coupling data.

```fortran
! Conceptual usage within a main program
program main_sim
  use couplings
  use fileio ! Assuming namdInfo is here
  implicit none

  type(overlap) :: all_couplings, selected_couplings
  type(namdInfo) :: sim_params

  ! ... Read sim_params from input file ...

  ! Allocate and load/calculate couplings
  call TDCoupIJ("path/to/rundirs", all_couplings, selected_couplings, sim_params)

  ! selected_couplings%Dij and selected_couplings%Eig now contain the necessary
  ! NACs and energies for the simulation within the band range specified by sim_params.

  ! ... proceed with simulation using selected_couplings ...

end program main_sim
```
If `COUPCAR` is not present and `sim_params%LCPTXT` is true, the program expects `EIGTXT` and either `NATXT` or (`NATXT_RE`, `NATXT_IM`) to exist. The original direct calculation from WAVECARs within `TDCoupIJ` is commented out, emphasizing reliance on pre-existing text files.

## Dependencies and Interactions
- **`prec`**: For numerical precision definitions.
- **`constants`**: For physical constants like `imgUnit`.
- **`lattice`**: For lattice-related information (not directly used in shown subroutines but likely by `wavecar`).
- **`wavecar` module**: Crucial for reading WAVECARs (`sysinfo`, `LOADWAVE`) when NACs are calculated directly (though this path is commented out in `TDCoupIJ`).
- **`fileio` module**: Likely for `namdInfo` type definition and possibly other I/O helpers.
- **Input files**: `COUPCAR` (binary), or `EIGTXT`, `NATXT`/`NATXT_RE`/`NATXT_IM` (text). WAVECAR files if direct calculation were active.
- **Output files**: Can create `COUPCAR` (if calculated and `CoupToFile` is called), or `EIGTXT`, `NATXT` (via `writeNaEig`), and `ACEIGTXT` (via `initspace`).
- The module expects WAVECARs to be from Gamma-point calculations (single k-point).
- The NAC calculation `CoupIJ` uses the formula `(<psi_i(t)|psi_j(t+dt)> - <psi_j(t)|psi_i(t+dt)>)`, which is proportional to the actual NAC; the division by `2*dt` and multiplication by constants might be handled elsewhere or implicitly.
```
