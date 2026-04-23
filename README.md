# TStep: Jacobian-Free Global Stability Analysis Engine

> **Operational Status**: Active Development. Currently validated against OpenFOAM, structurally designed for solver-agnostic deployment.

## 1. System Objective
TStep is a highly optimized, time-stepping framework designed for the global linear stability analysis of compressible flows. By circumventing the explicit construction of the massive flow Jacobian matrix, TStep leverages an Arnoldi iteration algorithm and Fréchet derivative stencils to extract dominant Ritz eigenpairs directly from the outputs of an external CFD solver.

## 2. Mathematical Framework
The continuous flow evolution is linearized around a mathematically steady base state, $q_0$. The Jacobian matrix $A = \frac{\partial F(q_0)}{\partial q}$ dictates the amplification or decay of linear disturbances. 

TStep computes the action of the exponential matrix operator on a Krylov vector $v$ using finite-difference Fréchet derivatives. For a first-order accurate approximation, the matrix-vector product is computed as:

$$A v \approx \frac{F(q_0 + \epsilon_0 v) - F(q_0)}{\epsilon_0}$$

**The Noise Floor Constraint:** The operational viability of the Arnoldi iteration depends entirely on the flow solver's discretization error ($\epsilon_s = ||F(q_0) - q_0||$). The user-defined disturbance magnitude ($\epsilon_0$) must be mathematically balanced above $\epsilon_s$ to ensure the extracted eigenmodes represent physical instabilities rather than numerical noise.

## 3. Architecture & Data Flow

```text
=============================================================================
                     TStep Architecture & Data Flow
=============================================================================

[1] INITIALIZATION & WARMUP (main.f90)
 |-- configurationRead()       -> Parses inputs/inputs.in
 |-- read_flowfield()          -> Loads raw SACRED data
 |-- run_simulation()          -> [CFD WARMUP] Digests raw data into steady state
 |-- promote_to_initial()      -> Extracts absolute base state (q0)

[2] THE DIAGNOSTIC GATE (The Noise Floor)
 |-- run_simulation()          -> Computes F(q0)
 |-- pack_state()              -> Measures absolute solver drift (eps_s)
 |   |-- [FATAL CHECK]         -> IF (eps_s > TOL) ABORT (Base flow drifting)
 |   |-- [OPTIMIZATION]        -> IF (frechet_order == 1) Cache F(q0)

[3] ARNOLDI ENGINE (Arnoldi.f90)
 |-- initial_disturbance()     -> Seeds random unit-norm v_1
 |-- DO j = 1, krylov_size
 |    |-- apply_linearized_operator() -> Approximates A*v via Frechet Stencil
 |    |    |-- [CFD RUN]       -> Computes F(q_i) for active stencil nodes
 |    |    |-- [CACHE HIT]     -> Bypasses CFD for alpha=0 if order=1
 |    |-- Modified Gram-Schmidt
 |    |-- Reorthogonalization pass
 |-- LAPACK DGEEV              -> Solves eigensystem of Hessenberg matrix H_m
 |-- Ritz Extraction           -> Lifts eigenvectors back to full fluid domain

[4] OUTPUT (write_eigendata.f90)
 |-- write_eigen_files()       -> Generates eigenvalues.dat, eigenvectors.dat
=============================================================================
```

## 4. Directory Structure & Module Hierarchy

```text
TStep/
├── bin/              # Build artifacts, Makefile, and executable
├── docs/             # Operational guides (OpenBLAS, Integration)
├── inputs/           # Configuration files (inputs.in)
├── output/           # Output CSV and eigenvalue files
├── src/              # Core Fortran 90 source code
└── tests/            # Unit tests
```

### Module Dependencies & Compilation Order
The core engine is built on a strict dependency tree to ensure memory safety and compilation integrity:

```text
main.f90
  ├─ accuracy
  ├─ error_handling
  ├─ setup
  ├─ call_CFD
  ├─ random_disturbance
  ├─ read_flow
  │   └─ OpenFOAM_IO
  ├─ write_flow
  │   └─ OpenFOAM_IO
  ├─ write_output
  ├─ frechet_stencil
  └─ Arnoldi
      └─ state_vector
```

### Detailed Module Specifications

**`main.f90` (Dispatcher)**
Orchestrates the entire pipeline. Enforces the strict noise-floor validation before authorizing entry into the Arnoldi loop.

**`Arnoldi.f90` (Mathematical Engine)**
Computes Ritz eigenvalues and eigenvectors.
* **Features:** Classical Gram-Schmidt with one mandatory reorthogonalization pass for numerical stability, handles complex conjugate pairs from DGEEV, flexible output sorting (magnitude, real, imaginary).
* **Error States:** `801` (Zero norm), `802/803` (Invalid dimensions/Krylov size), `804` (Allocation failure), `805` (LAPACK/DGEEV failure).

**`call_CFD.f90` (Execution Wrapper)**
Executes external command and monitors exit status.
* **Features:** Synchronous execution, blocks until CFD completion, verifies exit codes.
* **Error States:** `1` (Launch failed), `2` (Non-zero exit code from solver).

**`random_disturbance.f90` (State Seeding)**
Generates normalized, uniform random perturbation vectors scaled strictly by `eps_0` for sensitivity analysis and Krylov initialization.

**`OpenFOAM_IO.f90` (Native Format Interface)**
Low-level read/write handling for OpenFOAM scalar (`p`, `rho`, `T`) and vector (`U`, `C`) fields. Preserves native headers during write operations.

**`state_vector.f90` (Memory Translation)**
Flattens 3D spatial domains into contiguous 1D arrays for LAPACK matrix operations, and unpacks the modified state back into physical dimensions for the CFD solver.

**`error_handling.f90` (System Diagnostics)**
Centralized hierarchical logging. No module writes directly to `stdout` for errors; all pass through `log_error(code, context)`.

## 5. Configuration (`inputs.in`)

The engine is driven exclusively by namelists located in `inputs/inputs.in`.

```fortran
! General Settings (solver-agnostic)
&General
    flow_format = 'OpenFOAM',     
    output_file = '../output/flowfield.csv',
    COMMAND_RUN = 'cd /path/to/case && rhoCentralFoam',
/

! Arnoldi & Fréchet Parameters
&Arnoldi
    krylov_size        = 50,
    frechet_order      = 1,
    eps_0              = 1.0d-6,
    TTime              = 1.0,
    eigenvalue_sort_by = 'magnitude',
/

! OpenFOAM-Specific Pointers
&OpenFOAM
    N_HEADER_grid = 21,
    N_HEADER_var  = 21,
    file_grid_in  = '/path/to/openfoam/case/0/C',
    file_var_in   = '/path/to/openfoam/case/timestep/',
    file_var_out  = '/path/to/openfoam/case/output_timestep/',
/
```

### Parameter Dictionary
* `flow_format`: Target solver architecture.
* `COMMAND_RUN`: Exact shell string to execute the external CFD integration.
* `krylov_size`: Subspace dimension (number of eigenvalues requested).
* `frechet_order`: Truncation order (1 for forward difference, 2+ for central).
* `eps_0`: Disturbance magnitude (Must be $> \epsilon_s$).
* `TTime`: The physical integration time ($\tau$).
* `eigenvalue_sort_by`: Defines spectral output ordering (`magnitude`, `real`, `imaginary`).

## 6. Solver Integration & File Formats

### Native OpenFOAM Support
TStep expects standard OpenFOAM dictionary formatting for vectors and scalars:
```text
<header lines>
N
(
value1
value2
...
)
```

### Adding a New Solver
The architecture is designed to integrate with external solvers without modifying the core Arnoldi logic. To add a new solver (e.g., SU2, CGNS):
1.  Implement a specific I/O module matching the signature of `OpenFOAM_IO.f90`.
2.  Add required state variables to `variables.f90`.
3.  Create a specific namelist definition in `setup.f90`.
4.  Add the new configuration block to `inputs.in`.
5.  Route the new `flow_format` string within `read_flow.f90` and `write_flow.f90`.

## 7. Execution Protocol

**Dependencies:**
* Fortran 90/95 Compiler (`gfortran`, `ifort`).
* OpenBLAS / LAPACK libraries (`libblas-dev`, `liblapack-dev`).
* POSIX-compliant execution environment.

**Compilation:**
Via standard Make utility:
```bash
cd bin
make clean
make
```

Manual Compilation Mapping:
```bash
gfortran -c src/accuracy.f90 -o obj/accuracy.o -J mod/
gfortran -c src/error_handling.f90 -o obj/error_handling.o -J mod/
# ... [Compile dependencies in order] ...
gfortran obj/*.o -o bin/tstep -llapack -lblas
```

**Execution:**
```bash
./tstep
```

## 8. Hierarchical Error Handling

The system utilizes strict numeric exit codes. `ERROR_CODE = MODULE_ID × 100 + SPECIFIC_ERROR`.

| Range     | Originating Module        | Target Failure States                                  |
|-----------|---------------------------|--------------------------------------------------------|
| `1-99`    | `MAIN`                    | Base flow transient, configuration corrupt.            |
| `100-199` | `SETUP`                   | Bad namelist routing, missing param file.              |
| `200-299` | `OPENFOAM_IO (Scalars)`   | File unreadable, header mismatch, allocation fail.     |
| `300-399` | `OPENFOAM_IO (Vectors)`   | Premature EOF, vector parsing error.                   |
| `400-499` | `READ_FLOW`               | Field aggregation mismatch.                            |
| `500-599` | `WRITE_OUTPUT`            | Access denied to export directory.                     |
| `600-699` | `CALL_CFD`                | Solver binary not found, solver crashed.               |
| `700-799` | `RANDOM_DISTURBANCE`      | Zero norm collision, memory cap hit.                   |
| `800-899` | `ARNOLDI`                 | Krylov breakdown, LAPACK workspace failure.            |
| `900-999` | `WRITE_EIGENDATA`         | I/O block on spectral output files.                    |

## 9. Memory & Constraints

### Memory Management
* **Zero Leaks:** All dynamically allocated arrays (`V`, `H`, workspace matrices) are tracked, checked via `ALLOCATED()`, and cleared before programmatic termination via `cleanup_allocations()`.
* **In-Memory Load:** TStep currently holds the base physical state and the active Krylov vectors in RAM simultaneously. 

### Hard Numerical Constraints: The Noise Floor Gate
Before initiating the Arnoldi loop, `main.f90` forces an absolute evaluation of the CFD discretization error.

If $||F(q_0) - q_0|| / ||q_0|| > \text{EPS\_S\_TOL}$:
* **Result:** FATAL ABORT. 
* **Cause:** The base state is numerically drifting. It is an active transient.
* **Resolution:** You cannot linearize a moving target. Extend your external base-flow solver execution until the residuals drop below the defined precision floor.

## 10. Roadmap & Version History

### Current Phase: Foundation & Validation
* Implementation of Jacobian-free Fréchet evaluations.
* Stable LAPACK DGEEV integration with reorthogonalized Gram-Schmidt.
* Noise Floor Gate implementation for rigorous data validation.

### Version History
* **v3.0 (Current):** Full Arnoldi eigenvalue analysis. Native caching of base state evaluation. Implementation of the solver noise floor safety gate.
* **v2.1:** Centralized hierarchical error handling integration.
* **v2.0:** Zeus cluster integration preparation.
* **v1.1:** External command execution support.
* **v1.0:** Base IO architecture.
