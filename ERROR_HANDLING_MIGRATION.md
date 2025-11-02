# Error Handling System Migration Guide

## Overview

This document describes the new centralized error handling system and how to migrate existing code to use it.

## Error Code Structure

### Hierarchical System
```
ERROR_CODE = MODULE_ID × 100 + SPECIFIC_ERROR
```

### Module ID Ranges
| Range   | Module                  | Purpose                           |
|---------|-------------------------|-----------------------------------|
| 0-99    | Main Program            | Top-level execution errors        |
| 100-199 | Setup/Configuration     | Configuration file errors         |
| 200-299 | OpenFOAM I/O (Scalars)  | Scalar field reading errors       |
| 300-399 | OpenFOAM I/O (Vectors)  | Vector field reading errors       |
| 400-499 | Flow Reading            | High-level flow field errors      |
| 500-599 | Output Writing          | CSV output errors                 |
| 600-699 | External Commands       | Command execution errors          |
| 700-799 | Random Disturbance      | Perturbation generation errors    |

## Benefits

1. **Unique error codes**: No conflicts between modules
2. **Easy identification**: Error code reveals which module failed
3. **Scalability**: Each module has 100 error codes (plenty of room)
4. **Centralized documentation**: All errors defined in one place
5. **Consistent logging**: Uniform error message format
6. **Debugging**: Can quickly trace error source from exit code

## Migration Steps

### Step 1: Update Module Error Codes

**Old approach (setup.f90):**
```fortran
IF (status_id /= 0) THEN
    WRITE(*,*) 'ERROR: Could not open file'
    ierr = 1
    RETURN
END IF
```

**New approach:**
```fortran
USE error_handling

IF (status_id /= 0) THEN
    ierr = ERR_SETUP_FILE_OPEN
    CALL log_error(ERR_SETUP_FILE_OPEN, 'File: '//TRIM(filename))
    RETURN
END IF
```

### Step 2: Update Main Program

**Old:**
```fortran
IF (error_status /= 0) THEN
    WRITE(*,*) 'FATAL: Configuration read failed.'
    RETURN
END IF
```

**New:**
```fortran
USE error_handling

IF (error_status /= 0) THEN
    CALL log_error(ERR_MAIN_CONFIG)
    STOP ERR_MAIN_CONFIG
END IF
```

### Step 3: Update OpenFOAM_IO Module

Replace hardcoded constants with imports:

```fortran
MODULE OpenFOAM_IO
    USE accuracy
    USE variables
    USE setup, ONLY: get_unit
    USE error_handling  ! NEW
    
    IMPLICIT NONE
    
    ! Remove local error code definitions - now imported from error_handling
    
CONTAINS
    
    SUBROUTINE read_OF_scalars(...)
        ! ...
        IF (iostat_val /= 0) THEN
            ierr = ERR_SCALAR_OPEN
            CALL log_error(ERR_SCALAR_OPEN, 'File: '//TRIM(filename))
            RETURN
        END IF
        ! ...
    END SUBROUTINE
    
END MODULE OpenFOAM_IO
```

## Example Error Output

**Before (inconsistent):**
```
ERROR in read_OF_scalars: Could not open file: /path/to/file
```

**After (structured):**
```
=========================================
ERROR CODE: 201
MODULE:     OPENFOAM_IO_SCALAR
DESCRIPTION: Cannot open scalar file
DETAILS:    File: /path/to/file
=========================================
```

## Module-Specific Updates

### setup.f90
```fortran
USE error_handling

! Replace ierr = 1 with ierr = ERR_SETUP_FILE_OPEN
! Replace ierr = 2 with ierr = ERR_SETUP_NAMELIST_READ
```

### call_CFD.f90
```fortran
USE error_handling

! Replace ERROR_STATUS = 1 with ERROR_STATUS = ERR_CMD_LAUNCH_FAILED
! Replace ERROR_STATUS = 2 with ERROR_STATUS = ERR_CMD_NONZERO_EXIT
```

### random_disturbance.f90
```fortran
USE error_handling

! Replace ERROR_STATUS = 1 with ERROR_STATUS = ERR_DIST_INVALID_LENGTH
! Replace ERROR_STATUS = 2 with ERROR_STATUS = ERR_DIST_ALLOC
! Replace ERROR_STATUS = 3 with ERROR_STATUS = ERR_DIST_ZERO_NORM
```

### read_flow.f90
```fortran
USE error_handling

! Add specific error codes for each failure point
IF (error_status /= 0) THEN
    error_status = ERR_FLOW_PRESSURE
    CALL log_error(ERR_FLOW_PRESSURE, 'File: '//TRIM(file_var)//'p')
    RETURN
END IF
```

### write_output.f90
```fortran
USE error_handling

IF (iostat_val /= 0) THEN
    error_status = ERR_OUTPUT_FILE_OPEN
    CALL log_error(ERR_OUTPUT_FILE_OPEN, 'File: '//TRIM(output_file))
    RETURN
END IF
```

## Compilation Order

Since `error_handling` module depends only on `accuracy`, compile it early:

```bash
gfortran -c src/accuracy.f90 -o obj/accuracy.o -J mod/
gfortran -c src/error_handling.f90 -o obj/error_handling.o -J mod/
gfortran -c src/variables.f90 -o obj/variables.o -J mod/
# ... rest of modules
```

## Testing Strategy

1. **Phase 1**: Add error_handling module to build
2. **Phase 2**: Update one module at a time, starting with main
3. **Phase 3**: Run test cases and verify error messages
4. **Phase 4**: Update documentation

## Backward Compatibility

During migration, you can support both systems:
- New code uses named constants from `error_handling`
- Old code continues with numeric literals
- Gradually migrate module by module

## Future Enhancements

1. **Error logging to file**: Add optional file output
2. **Error recovery strategies**: Define recovery procedures per error type
3. **Error statistics**: Track which errors occur most frequently
4. **Internationalization**: Support multiple languages for error messages
5. **Stack trace**: Add call stack information for debugging

## Quick Reference Card

| Module | Old Code | New Code |
|--------|----------|----------|
| Main | `ierr = 1` | `ierr = ERR_MAIN_CONFIG` |
| Setup | `ierr = 1` | `ierr = ERR_SETUP_FILE_OPEN` |
| Scalar I/O | `ierr = 1` | `ierr = ERR_SCALAR_OPEN` |
| Vector I/O | `ierr = 11` | `ierr = ERR_VECTOR_OPEN` |
| Flow Read | Custom msg | `ierr = ERR_FLOW_PRESSURE` |
| Output | Custom msg | `ierr = ERR_OUTPUT_FILE_OPEN` |
| Ext. Cmd | `ERROR_STATUS = 1` | `ERROR_STATUS = ERR_CMD_LAUNCH_FAILED` |
| Disturbance | `ERROR_STATUS = 1` | `ERROR_STATUS = ERR_DIST_INVALID_LENGTH` |
