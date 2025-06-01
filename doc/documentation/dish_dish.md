# dish.f90 (from src/dish)

## Overview
The `dish` module implements the Decoherence Induced Surface Hopping algorithm. This is a method used in nonadiabatic molecular dynamics to simulate transitions between electronic states, incorporating quantum decoherence effects. The module manages the propagation of quantum amplitudes, calculation of decoherence times, stochastic selection of decoherence events, projection of the wavefunction, and averaging over multiple trajectories.

## Key Components
- **`module dish`**:
    - Uses modules `prec` (precision), `fileio` (for `namdInfo` type and possibly `TDKS` type), `hamil` (for Hamiltonian, though not directly called here, it's used by `TimeProp`), and `TimeProp` (for quantum state propagation).

- **`subroutine init_random_seed()`**:
    - Initializes the random number generator using the system clock to ensure different random sequences for different runs.

- **`subroutine calcdect(ks, inp, DECOTIME, COEFFISQ)`**:
    - Calculates decoherence times (`DECOTIME`) for each basis state.
    - `ks`: `TDKS` type, contains current wavefunction `psi_c`.
    - `inp`: `namdInfo` type, contains `DEPHMATR` (dephasing rate matrix).
    - `COEFFISQ`: Output, `|psi_c(j)|^2`.
    - `DECOTIME(i) = 1.0 / sum_j (inp%DEPHMATR(j, i) * |psi_c(j)|^2)`.

- **`subroutine whichtodec(tion, inp, DECOTIME, which, decmoment, shuffle)`**:
    - Decides which state (`which`) the system should decohere to.
    - `tion`: Current time step.
    - `DECOTIME`: Calculated decoherence times for each state.
    - `decmoment`: Array tracking accumulated time since last potential decoherence event for each state.
    - `shuffle`: An array of band indices, shuffled randomly at each call.
    - It iterates through the `shuffle`d states. If `DECOTIME(i) <= decmoment(i)` for a state `i`, that state `i` is chosen as `which`, and `decmoment(which)` is reset to 0. `decmoment` for other states is incremented by `inp%POTIM` in `mpirundish` after this call.

- **`subroutine projector(ks, inp, COEFFISQ, which, tion, cstat, iend, fgend)`**:
    - Performs the projection (collapse) of the wavefunction `ks%psi_c`.
    - `which`: The state to potentially collapse to, determined by `whichtodec`.
    - `cstat`: Input/Output, the current adiabatic state index.
    - `iend`: Target state for recombination counting.
    - `fgend`: Flag indicating if recombination to `iend` has occurred.
    - Calculates a Boltzmann-weighted probability: `popBoltz(which) = |psi_c(which)|^2 * exp(-dE/kbT)` if `dE > 0`. `dE` is the energy difference between the current state `cstat` and the target state `which`, adjusted based on `inp%LHOLE` (electron or hole dynamics).
    - A random number `r` is generated.
        - If `r <= popBoltz(which)`: The system collapses to state `which`. `ks%psi_c` becomes a vector with 1 at `which` and 0 elsewhere. `cstat` is updated to `which`. If `which` is the `iend` state and `fgend` was 0, recombination counters `ks%recom_pops` are updated.
        - Else (project out): `ks%psi_c(which)` is set to zero, and the wavefunction is renormalized. The logic for choosing a new `cstat` in this branch if `cstat == which` is commented out.

- **`subroutine runDISH(ks, inp)`**:
    - Main driver for the DISH simulation.
    - Initializes `ks%dish_pops` (time-dependent population of each state) and `ks%recom_pops` (recombination event populations).
    - Sets initial state `istat` and end state `iend` based on input parameters.
    - Calls `init_random_seed()`.
    - Loops `inp%NTRAJ` times (number of trajectories):
        - Calls `mpirundish` for each trajectory. This part is marked for potential OpenMP parallelization, but the comment indicates MPI is under development.
    - Averages `ks%dish_pops` over `inp%NTRAJ`. Sets initial population `ks%dish_pops(istat, 1) = 1.0_q`.

- **`subroutine mpirundish(ks, inp, istat, cstat, iend, fgend)`**:
    - Runs a single DISH trajectory. (The name "mpi" might be historical or for future MPI integration).
    - `istat`: Initial state for this trajectory.
    - `cstat`: Output, the final state of this trajectory.
    - Initializes `ks%psi_c` to be localized on `istat`. `cstat` is set to `istat`.
    - `decmoment` (accumulated time for decoherence check) is initialized to 0.
    - `shuffle` array is initialized.
    - Loops `tion` from 1 to `inp%RTIME - 1` (total simulation time steps):
        1.  `call PropagationT(ks, inp, tion)`: Propagates `ks%psi_c` for one time step.
        2.  `call calcdect(ks, inp, DECOTIME, COEFFISQ)`: Calculates decoherence times.
        3.  `call whichtodec(tion, inp, DECOTIME, which, decmoment_dummy, shuffle)`: Determines if a decoherence event should occur and to which state (`which`). (Note: `decmoment` is passed as a dummy argument here, the actual update happens below).
        4.  `decmoment = decmoment + inp%POTIM`: Increments accumulated time for all states.
        5.  If `which > 0` (a state was selected for decoherence):
            `call projector(ks, inp, COEFFISQ, which, tion, cstat, iend, fgend)`: Perform projection.
        6.  Updates population `ks%dish_pops(cstat, tion + 1)` by adding 1.0 for the current active state `cstat`.

- **`subroutine printDISH(ks, inp)`**:
    - Prints the simulation results.
    - Writes header information (simulation parameters) and time-dependent data to output files:
        - `SHPROP.NAMDTINI`: Main output. For each time step: time, average energy (`SUM(ks%eigKs*ks%dish_pops)`), and populations `ks%dish_pops(i,tion)`.
        - `RECOMB.NAMDTINI`: Recombination data. For each time step: time, average energy, and `ks%recom_pops(i,tion)`.
    - `NAMDTINI` is the initial time step number, used to create unique filenames if running multiple segments.

## Important Variables/Constants
- **`ks%psi_c`**: Complex array, wavefunction coefficients.
- **`ks%eigKs`**: Real array, Kohn-Sham eigenenergies at different times.
- **`ks%dish_pops`**: Real array (`NBASIS, RTIME`), stores the population of each state at each time step, averaged over trajectories.
- **`ks%recom_pops`**: Real array, similar to `dish_pops` but tracks populations related to recombination events.
- **`inp%DEPHMATR`**: Real matrix, dephasing rates between states, used in `calcdect`.
- **`inp%NTRAJ`**: Integer, number of trajectories to run.
- **`inp%RTIME`**: Integer, total number of time steps for the simulation.
- **`inp%POTIM`**: Real, ionic time step.
- **`inp%TEMP`**: Real, temperature (used for Boltzmann factor in `projector`).
- **`BOLKEV`**: Boltzmann constant in eV/K (likely from `prec` or `constants` module).
- **`cstat`**: Integer, current adiabatic electronic state the system is in.
- **`which`**: Integer, state selected for decoherence by `whichtodec`.
- **`decmoment`**: Real array, tracks time accumulation for decoherence decision.
- **`DECOTIME`**: Real array, calculated characteristic decoherence time for each state.

## Usage Examples
The main entry point for running a DISH simulation is the `runDISH` subroutine.

```fortran
program run_dish_simulation
  use dish
  use fileio ! For TDKS, namdInfo types
  implicit none

  type(TDKS) :: kohn_sham_system
  type(namdInfo) :: simulation_parameters

  ! ... Initialize kohn_sham_system (e.g., load energies, couplings) ...
  ! ... Initialize simulation_parameters (e.g., read from input file) ...

  call runDISH(kohn_sham_system, simulation_parameters)
  call printDISH(kohn_sham_system, simulation_parameters)

end program run_dish_simulation
```

## Dependencies and Interactions
- **`prec`**: For precision kinds.
- **`fileio`**: For `namdInfo` and `TDKS` derived types.
- **`hamil`**: Used by `TimeProp` which is called by `dish`.
- **`TimeProp`**: Contains `PropagationT` subroutine used for evolving the wavefunction.
- The algorithm relies on random numbers, initialized by `init_random_seed`.
- Output files `SHPROP.X` and `RECOMB.X` (where X is `inp%NAMDTINI`) are generated by `printDISH`.
```
