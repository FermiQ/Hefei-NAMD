# nac.py

## Overview
This script calculates Nonadiabatic Couplings (NACs) between electronic states using Kohn-Sham wavefunctions typically obtained from VASP (Vienna Ab initio Simulation Package) calculations. NACs are essential for simulating molecular dynamics involving transitions between electronic states (e.g., surface hopping). The script reads WAVECAR files from different time steps, performs phase correction and orthogonalization if specified, and computes the NACs using a finite difference approximation. It supports parallel computation for efficiency.

## Key Components
- **`orthogon(cic)`**:
    - Orthogonalizes a set of input coefficient vectors `cic` using Löwdin orthogonalization.
    - `cic`: A numpy array of complex coefficients representing wavefunctions.
    - Returns `cio`: An array of orthogonalized coefficients.
- **`init_wav_from_vaspwfc(wave0, gamma=True, icor=1, bmin=None, bmax=None, omin=None, omax=None, ikpt=1, ispin=1, dt=1.0)`**:
    - Initializes and saves a reference wavefunction (`cio.npy`) from a given WAVECAR file (`wave0`). This reference is used for phase correction in subsequent NAC calculations.
    - `wave0`: Path to the reference WAVECAR file.
    - `gamma`: Boolean, True if WAVECAR is from a Gamma-point calculation.
    - `icor`: Integer, specifies the correction scheme (1: phase correction based on previous step, 2: phase correction based on a fixed reference `cio.npy`, 12: orthogonalization then phase correction based on `cio.npy`).
    - `bmin`, `bmax`: Minimum and maximum band indices for NAC calculation.
    - `omin`, `omax`: Minimum and maximum band indices for orthogonalization/reference.
    - `ikpt`, `ispin`: K-point and spin index.
- **`nac_from_vaspwfc(waveA, waveB, gamma=True, icor=1, bmin=None, bmax=None, omin=None, omax=None, ikpt=1, ispin=1, dt=1.0)`**:
    - Calculates NACs between two sets of wavefunctions from `waveA` (at time t) and `waveB` (at time t+dt).
    - It uses the formula: `<psi_i(t)| d/dt |psi_j(t)> \approx (<psi_i(t)|psi_j(t+dt)> - <psi_i(t+dt)|psi_j(t)>) / (2dt)`.
    - `waveA`, `waveB`: Paths to WAVECAR files at time t and t+dt.
    - `dt`: Ionic time step in fs.
    - Other parameters are similar to `init_wav_from_vaspwfc`.
    - Returns `EnT` (average energies), `nacs` (NAC matrix), `pij` (overlap <i(t)|j(t+dt)>), `pji` (overlap <i(t+dt)|j(t)>).
- **`parallel_nac_calc(runDirs, nproc=None, gamma=False, icor=1, bmin=None, bmax=None, omin=None, omax=None, ikpt=1, ispin=1, dt=1.0)`**:
    - Performs NAC calculations in parallel for a series of WAVECAR files.
    - `runDirs`: A list of paths to directories, each containing a WAVECAR file for a time step.
    - `nproc`: Number of processors to use; defaults to all available CPUs.
    - Saves results (`eigXX.txt`, `nacXX_re.txt`, `nacXX_im.txt`) in the respective directories.

## Important Variables/Constants
- **`vaspwfc`**: A module (presumably custom or from a library) used to read data from VASP WAVECAR files.
- **`cio.npy`**: File storing the reference orthogonalized wavefunction coefficients. Created by `init_wav_from_vaspwfc` if `icor` is 2 or 12, and used by `nac_from_vaspwfc` for phase correction.
- **`cc1.npy`**: File storing phase correction factors. Used when `icor=1`.
- **Input WAVECAR files**: The script expects WAVECAR files to be located in sequentially named directories (e.g., `./0001/WAVECAR`, `./0002/WAVECAR`, etc.).
- **Output files**:
    - `eigXX.txt`: Contains band energies.
    - `nacXX_re.txt`: Contains the real part of NACs.
    - `nacXX_im.txt`: Contains the imaginary part of NACs.
    (XX in the filename corresponds to the `icor` value).

## Usage Examples
The script is executed from the command line. The main block (`if __name__ == '__main__':`) configures parameters like the range of time steps (`T_start`, `T_end`), reference time step (`T_ref`), number of processors (`nproc`), k-point, spin, band ranges, and correction scheme.

Example configuration in the script:
```python
if __name__ == '__main__':
    T_start = 1
    T_end   = 4000
    T_ref   = 1      # Reference time step for cio.npy generation (if icor=2 or 12)
    nproc   = 4      # Use 4 processors
    gamma   = False  # Not a Gamma-point calculation
    icor    = 12     # Orthogonalize and use reference cio.npy
    bmin    = 316    # Min band index for NAC
    bmax    = 350    # Max band index for NAC
    omin    = bmin   # Min band index for orthogonalization
    omax    = bmax   # Max band index for orthogonalization
    ikpt    = 1
    ispin   = 1

    WaveCars = ['./%04d/WAVECAR' % (ii + 1) for ii in range(T_start-1, T_end)]

    # Conditional removal of previous reference files
    if T_start == T_ref:
        if os.path.exists("cio.npy"):
            os.remove("cio.npy")
        if os.path.exists("cc1.npy") and icor == 1:
            os.remove("cc1.npy")

    # Initialize reference wavefunction if needed
    if icor == 2 or icor == 12:
        wave0 = './%04d/WAVECAR' % (T_ref)
        init_wav_from_vaspwfc(wave0, gamma, icor, bmin, bmax, omin, omax, ikpt, ispin)

    # Perform parallel NAC calculation
    parallel_nac_calc(WaveCars, nproc, gamma, icor, bmin, bmax, omin, omax, ikpt, ispin)
```

## Dependencies and Interactions
- **NumPy**: Extensively used for numerical operations, linear algebra (`np.linalg.eig`), and array manipulations.
- **SciPy**: No direct import visible, but `vaspwfc` or other dependencies might use it.
- **`vaspwfc` module**: Crucial for reading WAVECAR files. This is likely a specialized module for VASP output processing.
- **Python `os` and `sys` modules**: Used for file system operations (checking/removing files, path manipulation) and system interactions (exit).
- **Python `multiprocessing` module**: Used for parallelizing NAC calculations across multiple WAVECAR pairs.
- **Input Data**: Requires a series of WAVECAR files from a molecular dynamics simulation, typically generated by VASP.
- **Output Data**: Generates text files containing energies and NAC components for each time step analyzed. If `icor` involves a reference wavefunction (`cio.npy`) or phase factors (`cc1.npy`), these are also created/used.
```
