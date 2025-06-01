# TimeProp.f90 (from src/dish)

## Overview
This Fortran module, `TimeProp`, is responsible for the time evolution of quantum states in the 'dish' part of the simulation code. It integrates the time-dependent Schrödinger equation (TDSE). The primary method implemented is `PropagationT`, which uses a Trotter decomposition formula, noted as being robust and efficient, allowing for potentially larger electronic time steps. This scheme was proposed by Akimov & Prezhdo (JCTC 2014, 10, 789) and revised by Dr. Li Yunhai, specifically requiring real-valued off-diagonal Hamiltonian elements.

## Key Components
- **`module TimeProp`**:
    - Contains the subroutines for time propagation.
    - Uses modules `prec` (for precision kinds), `couplings` (likely for nonadiabatic couplings), and `hamil` (for Hamiltonian construction).

- **`subroutine PropagationT(ks, inp, tion)`**:
    - Implements the time propagation of the wavefunction coefficients (`ks%psi_c`).
    - **Arguments**:
        - `ks`: A derived type `TDKS` (likely Time-Dependent Kohn-Sham), passed as `intent(inout)`. It holds the wavefunction coefficients (`psi_c`) and Hamiltonian components (`ham_c`).
        - `inp`: A derived type `namdInfo`, passed as `intent(in)`. It contains simulation parameters like `POTIM` (ionic time step), `NELM` (number of electronic steps per ionic step), and `NBASIS` (number of basis states).
        - `tion`: An integer, `intent(in)`, representing the current ionic time step.
    - **Method**:
        1.  Calculates `edt = inp%POTIM / inp%NELM`, the electronic time step.
        2.  Iterates `tele` from 1 to `inp%NELM` (inner electronic steps).
        3.  Inside the loop:
            - Calls `make_hamil_rtime(tion, tele, ks, inp)` to construct the Hamiltonian for the current electronic step. (A commented-out alternative `make_hamil2` also exists).
            - Applies the Liouville-Trotter algorithm: `exp[i*(L_{ij}+L_i)*dt] \approx exp(i*L_{ij}*dt/2) * exp(i*L_i*dt) * exp(i*L_{ij}*dt/2)`.
                - **First L_ij part**: Iterates through pairs of basis states (jj, kk) and applies the off-diagonal part of the propagator:
                    ```fortran
                    phi = 0.5_q * edt * miuno * ks%ham_c(kk, jj) / hbar
                    cos_phi = cos(phi)
                    sin_phi = sin(phi)
                    cjj = ks%psi_c(jj)
                    ckk = ks%psi_c(kk)
                    ks%psi_c(jj) =  cos_phi * cjj + sin_phi * ckk
                    ks%psi_c(kk) = -sin_phi * cjj + cos_phi * ckk
                    ```
                    (Note: `miuno` is the imaginary unit `(0.0, -1.0)`). The comment "Changed the matrix structure for NAC. Now Ham_c(i,j) is stored as ks%ham_c(j,i)" indicates a specific storage convention for `ham_c`.
                - **L_i part**: Iterates through each basis state `jj` and applies the diagonal part of the propagator:
                    ```fortran
                    phi = edt * miuno * ks%ham_c(jj, jj) / hbar
                    ks%psi_c(jj) = ks%psi_c(jj) * exp(phi)
                    ```
                - **Second L_ij part**: Similar to the first, but iterates `jj` downwards and `kk` downwards from `NBASIS`.
        4.  After the electronic steps, it checks if the norm of `ks%psi_c` is conserved (between 0.99 and 1.01). If not, it prints an error and stops.

- **Commented-out `subroutine Propagation(ks, inp, tion)`**:
    - An older propagation scheme.
    - Appears to use a finite difference method (first or second order) for propagation.
    - Involves calls to `make_hamil` and `hamil_act`.

## Important Variables/Constants
- **`hbar`**: Planck's constant (value not shown in this file, likely imported from `prec` or globally defined).
- **`miuno`**: `complex(kind=q), parameter :: miuno = (0.0_q, -1.0_q)`, representing `-i`.
- **`ks%psi_c`**: Array holding the complex wavefunction coefficients. This is the primary quantity being propagated.
- **`ks%ham_c`**: Array holding the Hamiltonian matrix elements. It's constructed by `make_hamil_rtime` and used in `PropagationT`. The comment indicates `ks%ham_c(j,i)` stores H_ij.
- **`inp%POTIM`**: Ionic time step.
- **`inp%NELM`**: Number of electronic steps per ionic step.
- **`inp%NBASIS`**: Number of basis states/orbitals.
- **`edt`**: Electronic time step (`POTIM / NELM`).

## Usage Examples
This module is not called directly by the user but is a core part of the 'dish' simulation engine. The `PropagationT` subroutine is called within an outer loop over ionic time steps (`tion`).

```fortran
! Simplified conceptual usage within a larger simulation loop
module main_program
  use TimeProp
  use prec
  use types_module ! Assuming TDKS and namdInfo types are defined here
  implicit none

  type(TDKS) :: kohn_sham_state
  type(namdInfo) :: input_params
  integer :: current_ionic_step

  ! ... Initialize kohn_sham_state and input_params ...

  do current_ionic_step = 1, input_params%NAMDTIME - 1
    ! ... (Potentially update geometry, etc.) ...

    ! Propagate electronic wavefunction for one ionic step
    call PropagationT(kohn_sham_state, input_params, current_ionic_step)

    ! ... (Analyze results, surface hopping decisions, etc.) ...
  end do
end module main_program
```

## Dependencies and Interactions
- **`prec` module**: Provides precision kinds (e.g., `q` for real and complex numbers).
- **`couplings` module**: Likely provides routines or data related to nonadiabatic couplings, which would be part of the Hamiltonian construction.
- **`hamil` module**: Contains `make_hamil_rtime` (or `make_hamil2`) which constructs the Hamiltonian matrix used in the propagation.
- **`TDKS` derived type**: Defined elsewhere, encapsulates the quantum state (wavefunction coefficients, Hamiltonian).
- **`namdInfo` derived type**: Defined elsewhere, holds simulation setup parameters.
- The propagation scheme explicitly relies on the off-diagonal elements of `ks%ham_c` being real numbers.
- The norm of the wavefunction is checked after each full ionic step's worth of electronic propagations.

```
