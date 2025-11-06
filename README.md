# TStep - Time-Stepping Framework for CFD Solvers

> ⚠️ **Active Development**: This project is under continuous development. Features are being incrementally added and tested.

A Fortran-based time-stepping framework designed for general CFD solver integration. TStep is being developed as a flexible tool that can interface with various CFD solvers to perform time-accurate simulations, sensitivity analysis, and data processing.

## Development Status

**Current Phase**: Foundation and Testing with OpenFOAM  
**Target**: General-purpose time-stepping framework compatible with multiple CFD solvers

TStep is currently being developed and validated using OpenFOAM as the primary test case, but the architecture is designed to be **solver-agnostic**. The modular structure allows for easy adaptation to other CFD solvers (e.g., SU2, CFL3D, FUN3D, or custom solvers).

## Architecture Overview

TStep is built on a **modular, solver-agnostic architecture**:

1. **Solver Interface Layer**: Execute external CFD solvers and manage data exchange
2. **I/O Abstraction**: Read/write flow field data in various formats (currently OpenFOAM, extensible to others)
3. **Data Processing**: Manipulate flow fields (perturbations, filtering, etc.)
4. **Time-Stepping Core**: Control simulation advancement (under development)
5. **Configuration Management**: Flexible input system for different workflows

### Current Implementation: OpenFOAM Focus

While TStep is designed for general CFD solver compatibility, the current implementation focuses on OpenFOAM for validation and testing. The OpenFOAM-specific components are isolated in dedicated modules (`OpenFOAM_IO`) to facilitate future expansion to other solvers.

## Vision & Goals

**Ultimate Goal**: A general-purpose time-stepping framework that can:
- Interface with any CFD solver via external commands or API calls
- Support parallel execution and HPC environments

## Current Features

### ✅ Implemented
- **External command execution**: Run CFD simulations (OpenFOAM, or any solver) via shell commands
- **Random disturbance generation**: Create normalized, scaled random perturbation vectors for sensitivity analysis
- **OpenFOAM I/O**: Read and write scalar fields (pressure, density, temperature) and vector fields (velocity, coordinates)
- **Flow field perturbation**: Apply perturbations to flow fields and write back to OpenFOAM format
- **Data export**: Write processed data to CSV format
- **Flexible configuration**: Namelist-based input system
- **Robust error handling**: Centralized hierarchical error management system
- **Modular architecture**: Easy to extend for other solvers

### 🚧 In Development
- Solver abstraction layer for generic CFD integration
- Additional solver interfaces (SU2, custom solvers)
- Advanced perturbation methods
- Parallel execution support

### 📋 Planned
- Real-time monitoring and diagnostics

## Project Structure

```
TStep/
├── src/              # Source files
│   ├── main.f90              # Main program
│   ├── accuracy.f90          # Precision definitions
│   ├── error_handling.f90    # Centralized error code management
│   ├── variables.f90         # Global variables
│   ├── setup.f90             # Configuration reading
│   ├── call_CFD.f90          # External command execution
│   ├── random_disturbance.f90 # Random perturbation generation
│   ├── read_flow.f90         # Flow field reading module
│   ├── write_flow.f90        # Flow field writing module (OpenFOAM format)
│   ├── write_output.f90      # CSV output writing
│   └── OpenFOAM_IO.f90       # OpenFOAM file I/O routines (read & write)
├── inputs/           # Input configuration files
├── output/           # Output CSV files
├── bin/              # Compiled executables
├── mod/              # Compiled modules
└── obj/              # Object files
```

## Building

The project uses Fortran 90/95 with the following module dependencies:

```
main.f90
  ├─ accuracy
  ├─ error_handling
  │   └─ accuracy
  ├─ setup
  │   ├─ accuracy
  │   ├─ variables
  │   └─ error_handling
  ├─ call_CFD
  │   └─ error_handling
  ├─ random_disturbance
  │   ├─ accuracy
  │   ├─ variables
  │   ├─ setup
  │   └─ error_handling
  ├─ read_flow
  │   ├─ accuracy
  │   ├─ variables
  │   ├─ setup
  │   ├─ OpenFOAM_IO
  │   └─ error_handling
  ├─ write_flow
  │   ├─ accuracy
  │   ├─ variables
  │   ├─ setup
  │   ├─ OpenFOAM_IO
  │   └─ error_handling
  └─ write_output
      ├─ accuracy
      ├─ variables
      ├─ setup
      └─ error_handling
```

Compilation example (using gfortran):
```bash
gfortran -c src/accuracy.f90 -o obj/accuracy.o -J mod/
gfortran -c src/error_handling.f90 -o obj/error_handling.o -J mod/
gfortran -c src/variables.f90 -o obj/variables.o -J mod/
gfortran -c src/setup.f90 -o obj/setup.o -J mod/
gfortran -c src/call_CFD.f90 -o obj/call_CFD.o -J mod/
gfortran -c src/random_disturbance.f90 -o obj/random_disturbance.o -J mod/
gfortran -c src/OpenFOAM_IO.f90 -o obj/OpenFOAM_IO.o -J mod/
gfortran -c src/read_flow.f90 -o obj/read_flow.o -J mod/
gfortran -c src/write_flow.f90 -o obj/write_flow.o -J mod/
gfortran -c src/write_output.f90 -o obj/write_output.o -J mod/
gfortran -c src/main.f90 -o obj/main.o -J mod/
gfortran obj/*.o -o bin/TStep
```

## Solver Compatibility

### Current: OpenFOAM
Fully supported for development and testing. Reads and writes standard OpenFOAM file formats:
- Scalar fields: `p`, `rho`, `T`
- Vector fields: `U`, cell coordinates `C`

### Future: Generic Solver Support
The modular design allows for easy extension to other solvers:
- **File-based solvers**: Any solver that writes field data to files (Tecplot, CGNS, HDF5, etc.)
- **API-based solvers**: Direct in-memory coupling with solver libraries
- **Custom formats**: Add new I/O modules following the existing patterns

**Adding a new solver** requires:
1. Implement solver-specific I/O module (similar to `OpenFOAM_IO.f90`)
2. Add format option to configuration
3. Update `read_flow.f90` with new format case

## Configuration

The program reads configuration from `inputs/inputs.in` using a Fortran namelist format:

```fortran
&Setup
output_file         = '../output/flowfield.csv',
flow_format         = 'OpenFOAM',
COMMAND_RUN         = 'cd /path/to/case && rhoCentralFoam',
dist_mag            = 1.0d-6,
N_HEADER_grid       = 21,
N_HEADER_var        = 21,
file_grid_in        = '/path/to/openfoam/case/0/C',
file_var_in         = '/path/to/openfoam/case/timestep/',
file_var_out        = '/path/to/openfoam/case/output_timestep/',
/
```

### Configuration Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `output_file` | string | Path to output CSV file |
| `flow_format` | string | Flow solver format ('OpenFOAM', extensible to 'SU2', 'CGNS', etc.) |
| `COMMAND_RUN` | string | External command to execute CFD simulation |
| `dist_mag` | real | Magnitude for random disturbance vector scaling |
| `N_HEADER_grid` | integer | Number of header lines in grid coordinate file |
| `N_HEADER_var` | integer | Number of header lines in variable files |
| `file_grid_in` | string | Path to OpenFOAM cell center coordinate file (typically `C`) |
| `file_var_in` | string | Path prefix to OpenFOAM time directory containing input field variables |
| `file_var_out` | string | Path prefix to OpenFOAM time directory for writing output field variables |

## Input File Format

The program expects OpenFOAM files in the standard format:
```
<header lines>
N
(
value1
value2
...
valueN
)
```

For vector fields:
```
<header lines>
N
(
(x1 y1 z1)
(x2 y2 z2)
...
(xN yN zN)
)
```

## Output Format

The program generates a CSV file with the following columns:
```
X,Y,Z,rho,p,T,U,V,W
```

Where:
- `X, Y, Z`: Cell center coordinates
- `rho`: Density
- `p`: Pressure
- `T`: Temperature
- `U, V, W`: Velocity components

## Usage

1. Create or modify the configuration file `inputs/inputs.in`
2. Run the executable from the `bin/` directory:
   ```bash
   cd bin
   ./TStep
   ```
3. The output CSV file will be created at the location specified in the configuration

## Error Handling System

TStep implements a **centralized hierarchical error handling system** where each module owns a unique range of error codes.

### Error Code Structure

```
ERROR_CODE = MODULE_ID × 100 + SPECIFIC_ERROR
```

### Module Error Code Ranges

| Range   | Module                  | Example Errors                     |
|---------|-------------------------|---------------------------------|
| 1-99    | Main Program            | Configuration, flow read, output write failures |
| 100-199 | Setup/Configuration     | File open, namelist parse errors  |
| 200-299 | OpenFOAM I/O (Scalars)  | Scalar file read, allocation errors |
| 300-399 | OpenFOAM I/O (Vectors)  | Vector file read, format errors   |
| 400-499 | Flow Reading            | Pressure, density, temperature, velocity read failures |
| 500-599 | Output Writing          | Output file open, write errors    |
| 600-699 | External Commands       | Command launch, execution failures |
| 700-799 | Random Disturbance      | Invalid length, allocation, normalization errors |

### Benefits

1. **Unique identification**: Each error code immediately identifies its source module
2. **Structured logging**: Unified `log_error()` function provides consistent error output
3. **Scalability**: Each module has 100 error codes (room for growth)
4. **Debugging**: Exit codes directly trace back to failure source
5. **Maintainability**: All error codes centrally defined in `error_handling` module

### Example Error Output

```
=========================================
ERROR CODE: 201
MODULE:     OPENFOAM_IO_SCALAR
DESCRIPTION: Cannot open scalar file
DETAILS:    File: /path/to/file/p
=========================================
```

## Modules

### call_CFD
External command execution module for running CFD simulations:
- `run_simulation(COMMAND_STRING, ERROR_STATUS)`: Executes external command and monitors exit status

#### Features:
- Synchronous execution (waits for command completion)
- Command launch status checking (`CMDSTAT`)
- Exit code verification (`EXITSTAT`)
- Detailed error reporting

#### Error Codes:
- `0`: Success
- `1`: Command failed to launch (e.g., executable not found)
- `2`: Command executed but returned non-zero exit code

#### Usage:
The command specified in `COMMAND_RUN` configuration parameter is executed before reading flow field data. This allows the application to:
1. Run a CFD simulation
2. Wait for completion
3. Process the resulting flow field data

### random_disturbance
Random perturbation generation module for sensitivity analysis:
- `initial_disturbance(VECTOR_LENGTH, SCALING_CONSTANT, FINAL_VECTOR, ERROR_STATUS)`: Generates normalized random disturbance vector

#### Features:
- Generates uniform random values in [0, 1)
- L2-normalization (unit vector)
- Configurable magnitude scaling
- Robust allocation and validation

#### Algorithm:
1. Allocate vector of specified length
2. Fill with uniform random values using `RANDOM_NUMBER`
3. Calculate L2-norm using intrinsic `NORM2` function
4. Normalize to unit magnitude
5. Scale by `dist_mag` parameter

#### Error Codes:
- `0`: Success
- `1`: Invalid vector length (≤ 0)
- `2`: Memory allocation failed
- `3`: Zero or near-zero norm (cannot normalize)

### error_handling
Centralized error code management module providing:
- **Named error constants**: All error codes defined with descriptive names (e.g., `ERR_SCALAR_OPEN = 201`)
- **`log_error(error_code, [additional_info])`**: Unified error logging with structured output
- **`get_module_name(error_code)`**: Returns module name from error code
- **`get_error_description(error_code)`**: Returns human-readable error description

#### Key Features:
- Hierarchical error code system (MODULE_ID × 100 + ERROR_NUM)
- No error code conflicts between modules
- Optional additional context in error messages
- Easy to extend with new modules and error types

### accuracy
Defines precision for integer and real variables using ISO Fortran intrinsic types:
- `ik`: 32-bit integers (`int32`)
- `rk`: 64-bit real numbers (`real64`)

### variables
Global variables module containing:
- `N_HEADER_grid`: Number of header lines in grid file
- `N_HEADER_var`: Number of header lines in variable files
- `file_grid_in`: Input grid coordinate file path
- `file_var_in`: Input variable file path prefix
- `file_var_out`: Output variable file path prefix
- `output_file`: Output CSV file path
- `flow_format`: Flow solver format identifier
- `COMMAND_RUN`: External command string
- `dist_mag`: Disturbance magnitude scaling factor

### setup
Configuration reading module with subroutines:
- `configurationRead(ierr)`: Reads namelist from `inputs/inputs.in`
- `get_unit(u)`: Returns an available file unit number (10-99)

### OpenFOAM_IO
Low-level I/O routines for OpenFOAM file formats:
- `read_OF_scalars(filename, n_header_lines, data_vector, n_data_points, ierr)`: Reads scalar field data
- `read_OF_vectors(filename, n_header_lines, x_vector, y_vector, z_vector, n_data_points, ierr)`: Reads vector field data
- `write_OF_scalars(filename, n_header_lines, data_vector, n_data_points, ierr)`: Writes scalar field data preserving header/footer
- `write_OF_vectors(filename, n_header_lines, x_vector, y_vector, z_vector, n_data_points, ierr)`: Writes vector field data preserving header/footer

#### Error Codes
**Scalar Read Errors:**
- `ERR_SCALAR_OPEN (1)`: Cannot open file
- `ERR_SCALAR_HEADER (2)`: Error reading header
- `ERR_SCALAR_COUNT (3)`: Cannot read data count
- `ERR_SCALAR_INVALID_COUNT (4)`: Invalid data count (≤0)
- `ERR_SCALAR_ALLOC (5)`: Memory allocation failed
- `ERR_SCALAR_SKIP (6)`: Error skipping discard line
- `ERR_SCALAR_READ (7)`: Error reading data values

**Vector Read Errors:**
- `ERR_VECTOR_OPEN (11)`: Cannot open file
- `ERR_VECTOR_HEADER (12)`: Error reading header
- `ERR_VECTOR_COUNT (13)`: Cannot read data count
- `ERR_VECTOR_SKIP (14)`: Error skipping discard line
- `ERR_VECTOR_INVALID_COUNT (15)`: Invalid data count (≤0)
- `ERR_VECTOR_ALLOC (16)`: Memory allocation failed
- `ERR_VECTOR_EOF (17)`: Premature end of file
- `ERR_VECTOR_FORMAT (18)`: Format error parsing vector data

### read_flow
High-level flow field reading module:
- `read_flowfield(rho_in, p_in, T_in, U_in, V_in, W_in, Xgrid, Ygrid, Zgrid, data_count, error_status)`: Orchestrates reading of all flow field data

Reads the following files:
- `{file_var_in}p`: Pressure field
- `{file_var_in}rho`: Density field
- `{file_var_in}T`: Temperature field
- `{file_var_in}U`: Velocity field
- `{file_grid_in}`: Cell center coordinates

### write_flow
High-level flow field writing module:
- `write_flowfield(rho_out, p_out, T_out, U_out, V_out, W_out, data_count, error_status)`: Orchestrates writing of all flow field data back to OpenFOAM format

Writes the following files:
- `{file_var_out}p`: Pressure field
- `{file_var_out}rho`: Density field
- `{file_var_out}T`: Temperature field
- `{file_var_out}U`: Velocity field

### write_output
Output writing module:
- `write_flowfield_data(Xgrid, Ygrid, Zgrid, rho_in, p_in, T_in, U_in, V_in, W_in, data_count, error_status)`: Writes all data to CSV file

Output format: Scientific notation with 15 significant digits (`ES23.15E3`)

### main
Main program that:
1. Reads configuration from `inputs.in`
2. Executes external CFD command via `run_simulation()` if configured
3. Reads flow field data via `read_flowfield()`
4. Generates random disturbance vector via `initial_disturbance()`
5. Applies perturbations to flow fields and writes perturbed data to OpenFOAM format via `write_flowfield()`
6. Writes original (unperturbed) results to CSV via `write_flowfield_data()`
7. Performs cleanup of all allocated memory
8. Returns exit codes:
   - `0`: Success
   - `1`: Configuration read failure
   - `2`: External command failure
   - `3`: Flow field read failure
   - `4`: Disturbance generation failure
   - `5`: Output write failure

## Error Handling

The program uses a **centralized hierarchical error handling system**:

### Principles:
- **Centralized definitions**: All error codes defined in `error_handling` module
- **Hierarchical structure**: Error codes organized by module (100 codes per module)
- **Consistent status codes**: All subroutines use `error_status`/`ierr` output parameters
- **Structured logging**: Unified `log_error()` provides formatted error messages
- **Early returns**: Functions return immediately on error
- **Defensive programming**: All allocations checked; deallocations protected by `ALLOCATED()`
- **Comprehensive validation**: Input parameters validated before processing
- **Clean error paths**: Memory properly deallocated on all error conditions
- **Meaningful exit codes**: Main program exits with specific module-based error codes

### Error Code Ranges:
- **Main program (1-99)**: Configuration (1), external command (2), flow read (3), disturbance (4), output write (5)
- **Setup (100-199)**: File open (101), namelist read (102), invalid param (103)
- **OpenFOAM scalars (200-299)**: File operations (201-207)
- **OpenFOAM vectors (300-399)**: File operations (301-308)
- **Flow reading (400-499)**: Pressure (401), density (402), temperature (403), velocity (404), grid (405)
- **Output writing (500-599)**: File open (501), header write (502), data write (503)
- **External commands (600-699)**: Launch failed (601), non-zero exit (602)
- **Random disturbance (700-799)**: Invalid length (701), allocation (702), zero norm (703)

### Usage Example:

```fortran
USE error_handling

CALL some_operation(data, error_status)
IF (error_status /= 0) THEN
    CALL log_error(ERR_MODULE_SPECIFIC, 'Additional context here')
    STOP error_status
END IF
```

## Memory Management

The application follows best practices for memory management:

- **Safe allocation**: All `ALLOCATE` statements include `STAT=` checks
- **Safe deallocation**: All deallocations protected by `ALLOCATED()` checks
- **Centralized cleanup**: `cleanup_allocations()` subroutine handles all main arrays
- **Error path cleanup**: Memory deallocated on all error paths before early return
- **Pre-allocation cleanup**: Arrays deallocated before reallocation to prevent leaks
- **No memory leaks**: All allocated arrays properly tracked and freed

## Requirements

### Software
- Fortran 90/95 compiler (e.g., gfortran, ifort)
- CFD solver (currently tested with OpenFOAM)
- POSIX-compliant shell (for external command execution)

### System
- Sufficient memory for flow field data (all data loaded into memory)
- File system access for input/output operations

### Optional (for development)
- HPC environment (for cluster testing on Zeus)
- Version control (git recommended)

## Current Limitations

### Technical
- All data must fit in memory (no streaming)
- Single timestep processing (time-stepping under development)
- Serial execution only (parallel support planned)

### Solver-Specific
- **OpenFOAM**: Assumes standard file structure with parentheses
- Other solvers: Not yet implemented (but architecture supports them)

### Configuration
- Fixed configuration file location: `../inputs/inputs.in` relative to execution directory
- Limited output formats (CSV only, more formats planned)

## License

[Specify your license here]

## Contributing

This project is under active development. Contributions, suggestions, and feedback are welcome!

### Areas for Contribution
- Additional solver interfaces
- Time-stepping algorithms
- Parallel execution support
- Testing and validation
- Documentation improvements

## Contact & Support

For questions about the project or collaboration opportunities, please open an issue on GitHub.

## Development Roadmap

### Phase 1: Foundation (Current)
- ✅ Basic I/O with OpenFOAM
- ✅ External command execution
- ✅ Error handling system
- ✅ Random perturbation generation
- 🚧 Time-stepping framework

### Phase 2: Solver Abstraction
- 🔜 Generic solver interface
- 🔜 Additional solver support (SU2, custom formats)
- 🔜 Solver-agnostic data structures

### Phase 3: Advanced Features
- 📋 Parallel execution
- 📋 Adaptive time-stepping
- 📋 In-memory coupling

### Phase 4: Production
- 📋 HPC optimization
- 📋 Comprehensive testing suite
- 📋 Documentation and examples

**Legend**: ✅ Complete | 🚧 In Progress | 🔜 Next | 📋 Planned

## Repository Branches

- **main**: Stable release branch
- **external-commands**: Adds external command execution capability
- **zeus**: Active development branch (Zeus cluster integration and new features)

## Version History

- **v2.1 (zeus branch)**: 
  - Centralized hierarchical error handling system
  - Improved error diagnostics and logging
  - Enhanced code documentation
- **v2.0 (zeus branch)**: Zeus cluster integration
- **v1.1 (external-commands)**: Added external command execution support
- **v1.0 (Initial release)**: Basic OpenFOAM to CSV converter
