# TStep - OpenFOAM Flow Field Data Converter

A Fortran application for reading OpenFOAM flow field data and exporting it to CSV format for post-processing and analysis. The application can also execute external CFD simulations before processing the data.

## Overview

TStep reads OpenFOAM flow field data files (pressure, density, temperature, velocity, and grid coordinates) and exports them as a structured CSV file. The program is designed to work with OpenFOAM cell center data in the standard OpenFOAM file format.

## Features

- **External command execution**: Run CFD simulations (e.g., OpenFOAM solvers) directly from the application
- **Random disturbance generation**: Create normalized, scaled random perturbation vectors for sensitivity analysis
- Reads OpenFOAM scalar fields (pressure, density, temperature)
- Reads OpenFOAM vector fields (velocity components U, V, W)
- Reads grid cell center coordinates (X, Y, Z)
- Exports all data to a single CSV file
- Configurable via namelist input file
- Robust error handling and memory management with comprehensive error codes

## Project Structure

```
TStep/
├── src/              # Source files
│   ├── main.f90              # Main program
│   ├── accuracy.f90          # Precision definitions
│   ├── variables.f90         # Global variables
│   ├── setup.f90             # Configuration reading
│   ├── call_CFD.f90          # External command execution
│   ├── random_disturbance.f90 # Random perturbation generation
│   ├── read_flow.f90         # Flow field reading module
│   ├── write_output.f90      # CSV output writing
│   └── OpenFOAM_IO.f90       # OpenFOAM file I/O routines
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
  ├─ setup
  ├─ call_CFD
  ├─ random_disturbance
  │   ├─ accuracy
  │   ├─ variables
  │   └─ setup
  ├─ read_flow
  │   ├─ accuracy
  │   ├─ variables
  │   ├─ setup
  │   └─ OpenFOAM_IO
  └─ write_output
      ├─ accuracy
      ├─ variables
      └─ setup
```

Compilation example (using gfortran):
```bash
gfortran -c src/accuracy.f90 -o obj/accuracy.o -J mod/
gfortran -c src/variables.f90 -o obj/variables.o -J mod/
gfortran -c src/setup.f90 -o obj/setup.o -J mod/
gfortran -c src/call_CFD.f90 -o obj/call_CFD.o -J mod/
gfortran -c src/random_disturbance.f90 -o obj/random_disturbance.o -J mod/
gfortran -c src/OpenFOAM_IO.f90 -o obj/OpenFOAM_IO.o -J mod/
gfortran -c src/read_flow.f90 -o obj/read_flow.o -J mod/
gfortran -c src/write_output.f90 -o obj/write_output.o -J mod/
gfortran -c src/main.f90 -o obj/main.o -J mod/
gfortran obj/*.o -o bin/TStep
```

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
file_grid           = '/path/to/openfoam/case/0/C',
file_var            = '/path/to/openfoam/case/timestep/',
/
```

### Configuration Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `output_file` | string | Path to output CSV file |
| `flow_format` | string | Flow solver format (e.g., 'OpenFOAM') |
| `COMMAND_RUN` | string | External command to execute CFD simulation |
| `dist_mag` | real | Magnitude for random disturbance vector scaling |
| `N_HEADER_grid` | integer | Number of header lines in grid coordinate file |
| `N_HEADER_var` | integer | Number of header lines in variable files |
| `file_grid` | string | Path to OpenFOAM cell center coordinate file (typically `C`) |
| `file_var` | string | Path prefix to OpenFOAM time directory containing field variables |

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

### accuracy
Defines precision for integer and real variables using ISO Fortran intrinsic types:
- `ik`: 32-bit integers (`int32`)
- `rk`: 64-bit real numbers (`real64`)

### variables
Global variables module containing:
- `N_HEADER_grid`: Number of header lines in grid file
- `N_HEADER_var`: Number of header lines in variable files
- `file_grid`: Grid coordinate file path
- `file_var`: Variable file path prefix
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
- `{OF_file_var}p`: Pressure field
- `{OF_file_var}rho`: Density field
- `{OF_file_var}T`: Temperature field
- `{OF_file_var}U`: Velocity field
- `{OF_file_grid}`: Cell center coordinates

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
5. Writes results to CSV via `write_flowfield_data()`
6. Performs cleanup of all allocated memory
7. Returns exit codes:
   - `0`: Success
   - `1`: Flow field read failure
   - `2`: Output write failure
   - `3`: Disturbance generation failure

## Error Handling

The program uses a comprehensive error handling strategy:

### Principles:
- **Consistent status codes**: All subroutines use `error_status`/`ierr` output parameters
- **Early returns**: Functions return immediately on error with descriptive messages
- **Defensive allocation**: All allocations checked before use; deallocations check `ALLOCATED()` status
- **Comprehensive validation**: Input parameters validated before processing
- **Clean error paths**: Memory properly deallocated on all error conditions
- **Meaningful exit codes**: Main program exits with specific codes for each failure type

### Error Code Ranges:
- **Main program**: 0 (success), 1-3 (specific failures)
- **Configuration**: 0 (success), 1-2 (file/parse errors)
- **OpenFOAM scalars**: 0 (success), 1-7 (detailed I/O errors)
- **OpenFOAM vectors**: 0 (success), 11-18 (detailed I/O errors)
- **External commands**: 0 (success), 1-2 (launch/execution errors)
- **Random disturbance**: 0 (success), 1-3 (validation/allocation/normalization errors)

## Memory Management

The application follows best practices for memory management:

- **Safe allocation**: All `ALLOCATE` statements include `STAT=` checks
- **Safe deallocation**: All deallocations protected by `ALLOCATED()` checks
- **Centralized cleanup**: `cleanup_allocations()` subroutine handles all main arrays
- **Error path cleanup**: Memory deallocated on all error paths before early return
- **Pre-allocation cleanup**: Arrays deallocated before reallocation to prevent leaks
- **No memory leaks**: All allocated arrays properly tracked and freed

## Requirements

- Fortran 90/95 compiler (e.g., gfortran, ifort)
- OpenFOAM case data with standard file format
- Sufficient memory for flow field data (all data loaded into memory)

## Limitations

- All data must fit in memory (no streaming)
- Fixed CSV output format
- Assumes standard OpenFOAM file structure with parentheses
- Configuration file must be at `../inputs/inputs.in` relative to execution directory

## License

[Specify your license here]

## Author

[Specify author information here]

## Branches

- **main**: Stable release branch
- **external-commands**: Adds external command execution capability
- **zeus**: Development branch for Zeus cluster integration

## Version History

- v2.0 (zeus branch): Zeus cluster integration
- v1.1 (external-commands): Added external command execution support
- v1.0 (Initial release): Basic OpenFOAM to CSV converter
