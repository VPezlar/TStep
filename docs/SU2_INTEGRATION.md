# SU2 Integration in TStep

## Overview

TStep supports SU2 as a second CFD backend alongside OpenFOAM. The solver is
selected via `flow_format = 'SU2'` in the `&General` namelist of
`inputs/inputs.in`. All SU2-specific I/O lives in `src/SU2_IO.f90`; the
mathematical core (`Arnoldi.f90`, `state_vector.f90`, `frechet_stencil.f90`)
is completely untouched.

## Conservative ↔ Primitive Conversion

SU2 stores **conservative** variables; TStep's state vector is **primitive**.

### Conservative → Primitive (read path)

```
U = Momentum_x / Density
V = Momentum_y / Density
W = Momentum_z / Density
p = (gamma_gas - 1) * (Energy - 0.5 * Density * (U² + V² + W²))
T = p / (Density * R_gas)
rho = Density
```

### Primitive → Conservative (write path)

```
Density    = rho
Momentum_x = rho * U
Momentum_y = rho * V
Momentum_z = rho * W
Energy     = p / (gamma_gas - 1) + 0.5 * rho * (U² + V² + W²)
```

`gamma_gas` and `R_gas` are configured in the `&SU2` namelist.

## SU2 Restart File Format

SU2's `RESTART_ASCII` output:
- **Line 1**: comma-separated, double-quoted column names
  (e.g. `"PointID","x","y","z","Density","Momentum_x",...`)
- **Lines 2..N+1**: comma-separated numeric data rows
- **Trailing lines**: optional metadata block (non-numeric first token)

Column indices are discovered at runtime by parsing the header — never
hard-coded. Required columns: `x`, `y`, `z`, `Density`, `Momentum_x`,
`Momentum_y`, `Momentum_z`, `Energy`.

## Key Differences from OpenFOAM

| Aspect | OpenFOAM | SU2 |
|--------|----------|-----|
| Grid | Separate file (`baseflow_grid`) | Embedded in restart |
| Field files | One file per field (p, rho, T, U) | Single restart file |
| Endpoint | Time-stamped folder (e.g. `1.1/`) | Static file (`su2_solution_out`) |
| Variables | Primitive | Conservative (converted in SU2_IO) |

## Files Modified

- `src/variables.f90` — SU2 globals (`su2_config_file`, `su2_restart_in`,
  `su2_solution_out`, `gamma_gas`, `R_gas`)
- `src/setup.f90` — `&SU2` namelist read + SU2 branches in all five utility
  routines (`validate_paths`, `seed_stability_initial`,
  `promote_to_initial_state`, `find_endpoint_folder`, `clear_endpoint_folder`)
- `src/error_handling.f90` — SU2 I/O error codes (1200–1299)
- `src/read_flow.f90` — SU2 dispatch in `read_flowfield`
- `src/write_flow.f90` — SU2 dispatch in `write_flowfield`
- `src/main.f90` — endpoint `INQUIRE` checks gated on `flow_format`
- `bin/makefile` — `SU2_IO.o` added to the build

## New Files

- `src/SU2_IO.f90` — `read_SU2_solution` and `write_SU2_restart`

## Configuration

In `inputs/inputs.in`, set `flow_format = 'SU2'` in `&General` and provide
the `&SU2` namelist:

```
&SU2
    su2_config_file  = '/path/to/flow.cfg',
    su2_restart_in   = '/path/to/restart_flow_in.csv',
    su2_solution_out = '/path/to/restart_flow.csv',
    gamma_gas        = 1.4d0,
    R_gas            = 287.058d0,
    stability_dir    = '/path/to/Stab_TS/',
/
```

- `su2_restart_in`: the restart file TStep writes into (SU2 reads as input)
- `su2_solution_out`: the restart file SU2 writes (TStep reads back)
- `su2_config_file`: path to the SU2 `.cfg` configuration file
- `COMMAND_RUN` in `&General`: the shell command to invoke SU2
  (e.g. `'cd /path/to/case && SU2_CFD flow.cfg'`)
