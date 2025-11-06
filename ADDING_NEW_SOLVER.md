# Adding a New Solver to TStep

This guide demonstrates how to add support for a new CFD solver to TStep. We'll use SU2 as an example.

## Step 1: Add Solver-Specific Variables

Edit `src/variables.f90` to add variables specific to the new solver:

```fortran
MODULE variables
    USE accuracy
    IMPLICIT NONE

    ! General variables
    INTEGER(ik) :: N_HEADER_grid, N_HEADER_var
    CHARACTER(len=256) :: file_grid_in, file_var_in, file_var_out
    CHARACTER(len=256) :: output_file, flow_format, COMMAND_RUN
    REAL(rk) :: dist_mag

    ! OpenFOAM-specific variables
    ! (existing OpenFOAM variables)

    ! SU2-specific variables
    CHARACTER(len=256) :: su2_config_file, su2_mesh_file
    CHARACTER(len=256) :: su2_solution_file
    INTEGER(ik) :: su2_restart_iter

END MODULE variables
```

## Step 2: Create Configuration Namelist

Edit `src/setup.f90` to add a new namelist for the solver:

```fortran
SUBROUTINE configurationRead(ierr)
    ! ... existing code ...

    ! SU2-specific namelist
    NAMELIST / SU2 / su2_config_file, &
                     su2_mesh_file, &
                     su2_solution_file, &
                     su2_restart_iter

    ! ... existing file opening code ...

    ! Read General namelist
    READ(unit_num, nml=General, iostat=status_id, iomsg=message)
    ! ... error handling ...

    ! Read solver-specific namelist
    IF (TRIM(flow_format) == 'OpenFOAM') THEN
        READ(unit_num, nml=OpenFOAM, iostat=status_id, iomsg=message)
        ! ... error handling ...
    ELSE IF (TRIM(flow_format) == 'SU2') THEN
        READ(unit_num, nml=SU2, iostat=status_id, iomsg=message)
        IF (status_id /= 0) THEN
            ierr = ERR_SETUP_NAMELIST_READ
            CALL log_error(ERR_SETUP_NAMELIST_READ, 'SU2 namelist - '//TRIM(message))
            CLOSE(unit_num)
            RETURN
        END IF
    ELSE
        ! Unsupported format
        ierr = ERR_SETUP_INVALID_PARAM
        CALL log_error(ERR_SETUP_INVALID_PARAM, 'Unsupported flow_format: '//TRIM(flow_format))
        CLOSE(unit_num)
        RETURN
    END IF

    CLOSE(unit_num)
END SUBROUTINE configurationRead
```

## Step 3: Update Configuration File

Add a new section to `inputs/inputs.in`:

```fortran
! General Settings (solver-agnostic)
&General
    flow_format = 'SU2',
    output_file = '../output/flowfield.csv',
    dist_mag    = 1.0d-6,
    COMMAND_RUN = 'cd /path/to/case && SU2_CFD config.cfg',
/

! SU2-Specific Settings
&SU2
    su2_config_file   = '/path/to/su2/config.cfg',
    su2_mesh_file     = '/path/to/su2/mesh.su2',
    su2_solution_file = '/path/to/su2/solution.dat',
    su2_restart_iter  = 1000,
/
```

## Step 4: Create I/O Module

Create `src/SU2_IO.f90`:

```fortran
MODULE SU2_IO
    USE accuracy
    USE variables
    USE setup, ONLY: get_unit
    USE error_handling

    IMPLICIT NONE

CONTAINS

    SUBROUTINE read_SU2_scalars(filename, data_vector, n_data_points, ierr)
        CHARACTER(len=*), INTENT(in) :: filename
        REAL(rk), DIMENSION(:), ALLOCATABLE, INTENT(out) :: data_vector
        INTEGER(ik), INTENT(out) :: n_data_points
        INTEGER(ik), INTENT(out) :: ierr

        ! Implementation for reading SU2 scalar data
        ! ...
    END SUBROUTINE read_SU2_scalars

    SUBROUTINE read_SU2_vectors(filename, x_vector, y_vector, z_vector, n_data_points, ierr)
        CHARACTER(len=*), INTENT(in) :: filename
        REAL(rk), DIMENSION(:), ALLOCATABLE, INTENT(out) :: x_vector, y_vector, z_vector
        INTEGER(ik), INTENT(out) :: n_data_points
        INTEGER(ik), INTENT(out) :: ierr

        ! Implementation for reading SU2 vector data
        ! ...
    END SUBROUTINE read_SU2_vectors

END MODULE SU2_IO
```

## Step 5: Update Read/Write Modules

Edit `src/read_flow.f90`:

```fortran
SUBROUTINE read_flowfield(rho_in, p_in, T_in, U_in, V_in, W_in, &
                          Xgrid, Ygrid, Zgrid, data_count, error_status)
    ! ... existing declarations ...

    IF (flow_format == 'OpenFOAM') THEN
        ! Existing OpenFOAM read logic
        ! ...
    ELSE IF (flow_format == 'SU2') THEN
        ! SU2 read logic
        CALL read_SU2_scalars(TRIM(su2_solution_file)//'_p', p_in, data_count, error_status)
        ! ... more SU2 reads ...
    END IF

END SUBROUTINE read_flowfield
```

Similarly, update `src/write_flow.f90` for SU2 output.

## Step 6: Update Compilation

Add the new module to your compilation script:

```bash
gfortran -c src/SU2_IO.f90 -o obj/SU2_IO.o -J mod/
```

## Summary

The modular architecture makes adding new solvers straightforward:

1. ✅ Solver-specific variables → `variables.f90`
2. ✅ Configuration namelist → `setup.f90`
3. ✅ Configuration section → `inputs/inputs.in`
4. ✅ I/O routines → `SU2_IO.f90` (new file)
5. ✅ Flow read/write logic → `read_flow.f90` and `write_flow.f90`
6. ✅ Compilation → add to build process

The `&General` namelist remains unchanged, keeping the configuration clean and maintainable!
