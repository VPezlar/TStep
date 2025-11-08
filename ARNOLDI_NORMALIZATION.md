# Arnoldi Initial Vector Normalization

## Overview

The Arnoldi module now supports an optional flag to skip the initial vector normalization step. This is useful when the input vector has already been normalized to a specific magnitude by the `initial_disturbance` module.

## Important Mathematical Note

The Arnoldi iteration **requires** the initial vector v₁ to have unit norm (||v₁|| = 1) for the algorithm to work correctly. The Krylov basis must be orthonormal, and v₁ sets the scale for the entire basis.

**If v₁ is not normalized:**
- The Hessenberg matrix H will have incorrect scaling
- Eigenvalues (Ritz values) may be inaccurate
- Numerical stability may be compromised

## Usage

### Default Behavior (Automatic Normalization)

By default, the module **will normalize** the initial vector:

```fortran
CALL arnoldi_eigenvalues(v_init, m, frechet_order, eps_0, TTime, &
                        eigenvalues, eigenvectors, error_status)
```

This is the safest option and works for any input vector.

### Skip Normalization (Pre-normalized Vector)

If you have already normalized v_init to unit length, you can skip the normalization:

```fortran
CALL arnoldi_eigenvalues(v_init, m, frechet_order, eps_0, TTime, &
                        eigenvalues, eigenvectors, error_status, &
                        skip_normalization=.TRUE.)
```

**Requirements when using `skip_normalization=.TRUE.`:**
1. You **MUST** ensure ||v_init|| = 1 before calling
2. The module will print a warning and check the norm
3. If ||v_init|| differs from 1 by more than 1.0E-6, a warning will be printed

### Example with initial_disturbance Module

```fortran
PROGRAM example
    USE accuracy
    USE random_disturbance
    USE Arnoldi
    
    REAL(rk), DIMENSION(:), ALLOCATABLE :: v_init
    REAL(rk) :: magnitude
    COMPLEX(rk), DIMENSION(:), ALLOCATABLE :: evals, evecs
    INTEGER(ik) :: ierr
    INTEGER(ik) :: n
    
    n = 1000  ! System size
    magnitude = 1.0E-6_rk  ! Desired magnitude for CFD
    
    ! Generate normalized and scaled disturbance
    CALL initial_disturbance(n, magnitude, v_init, ierr)
    IF (ierr /= 0) STOP 'Disturbance generation failed'
    
    ! v_init is now normalized and scaled: ||v_init|| = magnitude
    ! For Arnoldi, we need ||v|| = 1, so normalize it
    v_init = v_init / NORM2(v_init)
    
    ! Now call Arnoldi with skip_normalization
    CALL arnoldi_eigenvalues(v_init, 10, 'First', 1.0E-6_rk, 1.0_rk, &
                            evals, evecs, ierr, &
                            skip_normalization=.TRUE.)
    
    IF (ierr /= 0) STOP 'Arnoldi failed'
    
    ! Use results...
    
END PROGRAM
```

## Workflow with CFD

When integrating with CFD solver:

1. **Generate disturbance** with `initial_disturbance`:
   - Vector is normalized and scaled to desired CFD magnitude
   - ||v|| = magnitude (e.g., 1.0E-6)

2. **Before calling Arnoldi**:
   - Normalize to unit length: `v_init = v_init / NORM2(v_init)`
   - Now ||v_init|| = 1

3. **Call Arnoldi** with `skip_normalization=.TRUE.`:
   - Avoids redundant normalization
   - Arnoldi uses v_init directly as v₁

4. **After Arnoldi**:
   - Eigenvectors can be scaled back to CFD magnitude if needed
   - Or use them at unit scale for subsequent iterations

## Warning Output

When using `skip_normalization=.TRUE.`, the module will print:

```
WARNING: Skipping initial vector normalization.
         User must ensure ||v_init|| = 1 for correct results.
         Current ||v_init|| = 1.000000E+00
```

If the norm significantly differs from 1:

```
WARNING: Skipping initial vector normalization.
         User must ensure ||v_init|| = 1 for correct results.
         Current ||v_init|| = 1.234567E-06
         WARNING: ||v_init|| significantly differs from 1!
                  Results may be inaccurate.
```

## API Update

### Old Signature
```fortran
SUBROUTINE arnoldi_eigenvalues(v_init, m, frechet_order, eps_0, TTime, &
                               eigenvalues, eigenvectors, ERROR_STATUS)
```

### New Signature
```fortran
SUBROUTINE arnoldi_eigenvalues(v_init, m, frechet_order, eps_0, TTime, &
                               eigenvalues, eigenvectors, ERROR_STATUS, &
                               skip_normalization)
! skip_normalization - LOGICAL, INTENT(IN), OPTIONAL
!                      Default: .FALSE. (normalize v_init)
!                      If .TRUE.: skip normalization (user must ensure ||v_init|| = 1)
```

## Recommendations

1. **For testing/debugging**: Use default normalization (don't pass the flag)

2. **For production CFD runs**: 
   - Pre-normalize your vector
   - Use `skip_normalization=.TRUE.` to avoid redundant computation
   - Always verify ||v|| = 1 before calling

3. **When in doubt**: Let the module normalize (default behavior)

## Performance Note

The normalization step is very cheap computationally:
- One `NORM2()` call: O(n)
- One vector division: O(n)

Skipping it provides minimal performance benefit, but is useful for:
- Clarity in code workflow
- Ensuring specific scaling is preserved
- Avoiding floating-point rounding in double normalization
