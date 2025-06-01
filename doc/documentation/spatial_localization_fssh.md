# spatial_localization_fssh.py

## Overview
This script analyzes the spatial localization of Kohn-Sham (KS) orbitals and its interplay with nonadiabatic molecular dynamics, specifically using results from VASP (PROCAR files) and FSSH (Fewest Switches Surface Hopping) simulations. It calculates the contribution of specified atoms to each KS orbital, visualizes energy bands colored by this localization, and further processes FSSH state populations to understand the time evolution of charge localization on selected atomic groups.

## Key Components
- **`WeightFromPro(infile='PROCAR', whichAtom=None, spd=None)`**:
    - Extracts KS orbital energies and the sum of orbital weights (projections) onto specified atoms from a VASP PROCAR file.
    - `infile`: Path to the PROCAR file.
    - `whichAtom`: A list of atom indices (1-based in input, converted to 0-based internally) for which to sum the weights. If `None`, sums weights over all atoms.
    - `spd`: A list of integers specifying which s, p, d orbital contributions to sum (e.g., `[0,1,2,3]` for s+p). If `None`, sums total contribution.
    - Returns `Energies` (nspin, nkpts, nbands) and `Weights` (nspin, nkpts, nbands) which is the sum of projections onto `whichAtom`.
- **`parallel_wht(runDirs, whichAtoms, nproc=None)`**:
    - Reads PROCAR files from multiple directories (`runDirs`) in parallel to extract energies and atomic weights using `WeightFromPro`.
    - `runDirs`: List of simulation directory paths.
    - `whichAtoms`: List of atom indices passed to `WeightFromPro`.
    - `nproc`: Number of processors for parallel execution.
    - Returns aggregated `enr` (energies) and `wht` (weights) as NumPy arrays.

## Important Variables/Constants
- **PROCAR Processing & Initial Plotting:**
    - `nsw`: Number of simulation steps for PROCAR analysis.
    - `dt`: Time step (fs).
    - `prefix`: Directory prefix for PROCAR files (e.g., `../../waveprod/run`).
    - `runDirs`: List of paths to PROCAR containing directories.
    - `whichS`, `whichK`: Specific spin and k-point index to analyze.
    - `whichA`: NumPy array of atom indices (0-based) defining a group of interest (e.g., a molecule or fragment).
    - `Alabel`, `Blabel`: Labels for groups of atoms if comparing localization (e.g., for a donor/acceptor system, though the example focuses on `Alabel`).
    - `all_wht.npy`, `all_en.npy`: Files to save/load calculated weights and energies from PROCARs.
    - `ksen_wht.png`: Output plot showing energy bands colored by localization on `whichA`.
- **FSSH Data Processing & Second Plotting:**
    - `bmin`, `bmax`: Min and max band indices relevant to the FSSH dynamics.
    - `namdTime`: Duration of the FSSH simulation in time steps.
    - `potim`: Time step in FSSH simulation (fs).
    - `inpFiles`: Glob pattern for FSSH output files (e.g., `../../317/SHPROP.*`), which contain state coefficients.
    - `weight_sh.dat`: Output file storing the time-evolved weighted sum of FSSH state populations by their localization (`rho`).
    - `na.npy`, `ad.npy`, `et.npy`: Files to save/load calculated time derivatives and expectation values related to localization and FSSH dynamics.
        - `et`: Expectation value of localization <P(t) * L> where P is population and L is localization.
        - `na`: Nonadiabatic term related to change in population and localization.
        - `ad`: Adiabatic term related to change in localization for given populations.
    - `kspat.png`: Output plot showing `et`, cumulative `na`, and cumulative `ad` over time.

## Usage Examples
The script is run from the command line:
```bash
python spatial_localization_fssh.py
```
**Setup:**
1.  **PROCARs:** Expects PROCAR files in directories like `../../waveprod/run/0001/`, `../../waveprod/run/0002/`, ...
2.  **FSSH outputs:** Expects FSSH `SHPROP.*` files in a directory like `../../317/`.
3.  Modify script parameters (`nsw`, `prefix`, `whichA`, `bmin`, `bmax`, `inpFiles`, etc.) as needed for the specific dataset.

**Execution Flow:**
1.  **Localization from PROCARs:**
    - If `all_wht.npy` and `all_en.npy` don't exist, it calls `parallel_wht` to process PROCARs, calculating energies (`Enr`) and summed weights (`Wht`) for atoms in `whichA`. Results are saved.
    - Generates `ksen_wht.png`: A plot of energy bands (`Enr`) vs. time (step number). Bands are either sized (Method 1, commented out) or colored (Method 2, active) by their localization (`Wht`) on the atoms defined in `whichA`. Method 3 (LineCollection) is also an option but commented out.
2.  **Localization Dynamics with FSSH:**
    - If `SHPROP.*` files are found and `weight_sh.dat` doesn't exist:
        - Loads FSSH state populations (`psi_a`).
        - Takes the localization values (`Wht`) for the relevant band range (`bmin:bmax`).
        - Calculates `rho`: average localization of the FSSH active state(s) over time.
        - Calculates `dcdt` (time derivative of localization of KS states) and `dpdt` (time derivative of FSSH populations).
        - Calculates `NA` (nonadiabatic contribution to localization change), `AD` (adiabatic contribution), and `ET` (total localization expectation value). These are saved to `.npy` files.
    - If `et.npy`, `na.npy`, `ad.npy` exist (or were just computed):
        - Generates `kspat.png`: A plot showing `ET`, `cumsum(NA)`, and `cumsum(AD)` vs. time, illustrating different contributions to the evolution of charge localization.

## Dependencies and Interactions
- **NumPy**: Essential for all numerical operations, array handling, and loading/saving `.npy` files.
- **`os` module**: For file existence checks (`os.path.isfile`).
- **`re` module**: For parsing PROCAR files.
- **`glob` module**: For finding FSSH output files.
- **Matplotlib**: Used for generating all plots (`ksen_wht.png`, `kspat.png`).
    - `pyplot`: Core plotting.
    - `LineCollection`: For advanced line plotting (Method 3).
    - `mpl_toolkits.axes_grid1.make_axes_locatable`: For colorbar placement.
- **`multiprocessing` module**: For parallelizing PROCAR file processing.
- **Input Files**:
    - VASP PROCAR files.
    - FSSH `SHPROP.*` output files.
    - Potentially pre-computed `.npy` files (`all_wht.npy`, `all_en.npy`, `na.npy`, `ad.npy`, `et.npy`) and `weight_sh.dat`.
- **Output Files**:
    - `ksen_wht.png`: Plot of energy bands vs. time, colored by localization.
    - `kspat.png`: Plot of localization dynamics components.
    - `all_wht.npy`, `all_en.npy`: Cached PROCAR data.
    - `weight_sh.dat`, `na.npy`, `ad.npy`, `et.npy`: Cached FSSH processed data.
```
