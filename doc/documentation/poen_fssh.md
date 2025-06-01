# poen_fssh.py

## Overview
This script is designed to process and visualize data from molecular dynamics simulations, specifically focusing on results from VASP (PROCAR files) and FSSH (Fewest Switches Surface Hopping) calculations. It extracts band energies, calculates average energies over multiple runs, and correlates them with state population dynamics from FSSH simulations. The final output is a plot showing energy levels over time, colored by state population, along with the average energy of a particular state (e.g., average hole energy).

## Key Components
- **`EnergyFromPro(infile='PROCAR')`**:
    - Extracts band energies from a VASP PROCAR file.
    - `infile`: Path to the PROCAR file.
    - Returns a NumPy array `energies` with dimensions (nspin, nkpts, nbands).
- **`parallel_energy(runDirs, nproc=None)`**:
    - Reads PROCAR files from multiple directories (`runDirs`) in parallel to extract energies.
    - `runDirs`: A list of directory paths, each expected to contain a PROCAR file.
    - `nproc`: Number of processors to use for parallel processing. Defaults to the total number of CPU cores.
    - Returns a NumPy array containing energies aggregated from all PROCAR files.

## Important Variables/Constants
- **Input Data (Main Script Block):**
    - **`nsw`**: Number of simulation steps (e.g., 7000). Used to define `runDirs`.
    - **`prefix`**: Prefix for the directories containing PROCAR files (e.g., 'NAMD/run/').
    - **`runDirs`**: List of paths to simulation directories (e.g., `['NAMD/run/0001', 'NAMD/run/0002', ...]`).
    - **`all_en.npy`**: A NumPy binary file. If it exists, energies are loaded from it. Otherwise, energies are calculated using `parallel_energy` and then saved to this file.
    - **`bmin`, `bmax`**: Minimum and maximum band indices to consider from the energy data (e.g., 152, 236).
    - **`namdTime`**: Duration of the NAMD (FSSH) simulation in time steps (e.g., 6000).
    - **`potim`**: Time step duration (e.g., 1.0 fs), though not explicitly used in calculations shown, it's a common simulation parameter.
    - **`inpFiles`**: A list of FSSH output files, typically matching a pattern like 'NAMD/SHPROP.*'. These files contain time, current state energy, and state coefficients (populations).
- **Processed Data:**
    - **`po.npy`**: A NumPy binary file. If it exists, processed time, energy, and coefficient data are loaded from it. Otherwise, it's created by processing `inpFiles` and `all_en.npy`. It stores:
        - `Time`: Time array for each band.
        - `En_t`: Energy of each band over time, averaged over `iniTimes` and referenced to `EVBM`.
        - `Ci_t`: Coefficients (populations) of each band over time, averaged over `iniTimes`.
    - **`average_energy.dat`**: Text file storing the time and average energy of the current state from FSSH trajectories.
    - **`EVBM`**: Average energy of the valence band maximum (or highest considered band), used as an energy reference.
- **Plotting:**
    - **`kpoen.png`**: The output image file name for the generated plot.

## Usage Examples
The script is executed directly from the command line:
```bash
python poen_fssh.py
```
It expects:
1.  PROCAR files to be present in directories specified by `prefix` and `nsw` (e.g., `NAMD/run/0001/PROCAR`, ...). Alternatively, a pre-computed `all_en.npy` file.
2.  FSSH output files (e.g., `NAMD/SHPROP.1`, `NAMD/SHPROP.2`, ...) in the directory specified by `inpFiles` pattern.

The script will then:
- Load or calculate energies from PROCARs.
- Load or process FSSH data.
- Generate a plot `kpoen.png` visualizing energy levels colored by their populations over time.
- Save processed data to `po.npy` and `average_energy.dat`.

## Dependencies and Interactions
- **NumPy**: Heavily used for array manipulation, loading/saving binary data (`.npy`), and numerical calculations.
- **`os` module**: For file system interactions like checking if a file exists (`os.path.isfile`).
- **`re` module**: For regular expression-based parsing of PROCAR files.
- **`glob` module**: For finding files matching a pattern (e.g., `NAMD/SHPROP.*`).
- **`multiprocessing` module**: For parallelizing the reading of PROCAR files.
- **Matplotlib**: Used for generating the final plot.
    - `mpl.use('agg')`: Sets the backend to 'agg' for generating figures without a display server.
    - `pyplot`: For plotting functions.
    - `mpl_toolkits.axes_grid1.make_axes_locatable`: For managing colorbar layout.
- **Input Files**:
    - PROCAR files from VASP.
    - FSSH output files (e.g., `SHPROP.*`).
    - Potentially `all_en.npy` and `po.npy` if previously generated.
- **Output Files**:
    - `kpoen.png` (plot).
    - `all_en.npy` (aggregated energies from PROCARs).
    - `po.npy` (processed time, energy, and coefficient data for plotting).
    - `average_energy.dat` (average energy of the FSSH active state over time).
```
