# SYSTEM DIRECTIVE: TStep Engine - SU2 Abstraction Implementation

## MISSION OBJECTIVE
You are an expert Fortran 90 software architect. Your objective is to expand the `TStep` global stability analysis framework. Currently, the I/O and state-management layers are tightly coupled to OpenFOAM idiosyncrasies. You will implement support for the `SU2` CFD solver. 

## STRICT ARCHITECTURAL CONSTRAINTS
1. **Mathematical Purity:** You will NOT modify `Arnoldi.f90` or `state_vector.f90`. The mathematical core operates strictly on 1D vectors of primitive variables ($\rho, p, T, U, V, W$). It must remain completely blind to the external solver.
2. **The Border Wall:** SU2 operates on conservative variables ($\rho, \rho U, \rho V, \rho W, \rho E$). The translation between SU2's conservative variables and TStep's primitive variables MUST happen exclusively inside a newly created `SU2_IO.f90` module.
3. **No Dynamic Folders:** OpenFOAM creates dynamic time-step folders (e.g., `1.09999/`). SU2 overwrites static files (e.g., `restart_flow.dat`). You must abstract the directory hunting logic.

Execute the following four phases precisely.

---

### PHASE 1: Configuration Expansion
Modify the configuration definitions to accept SU2 parameters.

**Target 1: `variables.f90`**
* Add the following string variables to hold the SU2 paths: `su2_config_file`, `su2_restart_in`, `su2_solution_out`.
* Add thermodynamic constants required for state translation: `gamma_gas` (default 1.4) and `R_gas` (default 287.05).

**Target 2: `setup.f90`**
* Inside `configurationRead()`, define a new namelist `&SU2` containing the variables above.
* Add an `ELSE IF (TRIM(flow_format) == 'SU2')` block to read the `&SU2` namelist from the `inputs.in` file.

---

### PHASE 2: Structural Abstraction (The Gates)
Bypass the OpenFOAM directory hunting for static solvers.

**Target 1: `setup.f90` -> `find_endpoint_folder`**
* Gate the existing POSIX shell pipeline (`ls | awk | sort | tail`) inside an `IF (TRIM(flow_format) == 'OpenFOAM')` block.
* Add an `ELSE IF (TRIM(flow_format) == 'SU2')` block. For SU2, there is no folder to hunt. Assign the path string directly to the predefined `su2_solution_out` directory.

**Target 2: `setup.f90` -> `clear_endpoint_folder`**
* Gate the existing `rm -rf` logic inside the OpenFOAM block.
* For SU2, execute a command to delete the existing static endpoint file (e.g., `rm -f <su2_solution_out>/restart_flow.dat`) to guarantee a clean slate before the solver runs.

---

### PHASE 3: The Translation Border (SU2_IO.f90)
Create a new module named `SU2_IO.f90` from scratch. This module serves as the firewall. It reads conservative data, translates it to primitive data for TStep, and translates primitive data back to conservative data for the solver.

**Target 1: `read_SU2_solution` subroutine**
* **Input:** File path.
* **Output:** `rho`, `p`, `T`, `U`, `V`, `W` arrays (primitive variables).
* **Logic:**
  1. Parse the SU2 output file (CSV or DAT).
  2. Read the conservative columns: Density ($\rho$), Momentum ($\rho U, \rho V, \rho W$), and Energy ($\rho E$).
  3. Execute the translation to primitive:
     * Velocity: $U = (\rho U)/\rho$ (same for V, W)
     * Pressure: $p = (\gamma_{gas} - 1) \cdot [\rho E - 0.5 \cdot \rho(U^2 + V^2 + W^2)]$
     * Temperature: $T = p / (\rho \cdot R_{gas})$

**Target 2: `write_SU2_restart` subroutine**
* **Input:** File path, primitive arrays (`rho`, `p`, `T`, `U`, `V`, `W`).
* **Logic:**
  1. Translate primitive back to conservative:
     * Momentum: $(\rho U) = \rho \cdot U$
     * Energy: $\rho E = \frac{p}{\gamma_{gas} - 1} + 0.5 \cdot \rho(U^2 + V^2 + W^2)$
  2. Write the arrays out in the exact formatting required by an SU2 restart file. (Preserve the native SU2 headers).

---

### PHASE 4: The Dispatchers
Route the high-level commands through your new abstraction layer.

**Target 1: `read_flow.f90`**
* Inside `read_flowfield`, add an `ELSE IF (TRIM(flow_format) == 'SU2')` block.
* Call `read_SU2_solution(...)` from `SU2_IO` to populate the `rho_in, p_in, T_in, U_in, V_in, W_in` arrays.

**Target 2: `write_flow.f90`**
* Inside `write_flowfield`, add the SU2 gate.
* Call `write_SU2_restart(...)` to dump the perturbed state back to the disk.

## EXECUTION ORDER
Provide the complete Fortran 90 code modifications for **Phase 1 and Phase 2** first. Wait for my confirmation before proceeding to write the `SU2_IO.f90` module in Phase 3.