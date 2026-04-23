# Arnoldi Module - Eigenvalue Computation API

## Overview

The Arnoldi module has been completely rewritten to compute Ritz eigenvalues and eigenvectors. The module uses LAPACK for eigenvalue decomposition and integrates with the TStep configuration system.

## Module API

### `arnoldi_eigenvalues(v_init, m, frechet_order, eps_0, TTime, eigenvalues, eigenvectors, ERROR_STATUS)`

Computes Ritz eigenvalues and eigenvectors of a matrix operator using Arnoldi iteration.

**Input Parameters:**
- `v_init` (REAL(rk), DIMENSION(:), INTENT(IN))
  - Initial disturbance vector
  - Size determines system dimension n
  - Must be non-zero

- `m` (INTEGER(ik), INTENT(IN))
  - Krylov subspace size
  - Must satisfy: 0 < m ≤ n
  - Read from `inputs.in` via `krylov_size` parameter

- `frechet_order` (CHARACTER(len=*), INTENT(IN))
  - Order of Frechet derivative ('First', 'Second', etc.)
  - Currently unused (reserved for future)
  - Read from `inputs.in`

- `eps_0` (REAL(rk), INTENT(IN))
  - Epsilon parameter for perturbations
  - Currently unused (reserved for future)
  - Read from `inputs.in`

- `TTime` (REAL(rk), INTENT(IN))
  - Time parameter
  - Currently unused (reserved for future)
  - Read from `inputs.in`

**Output Parameters:**
- `eigenvalues` (COMPLEX(rk), DIMENSION(:), ALLOCATABLE, INTENT(OUT))
  - Ritz eigenvalues (size m)
  - Sorted by **descending imaginary part**: Im(λ₁) ≥ Im(λ₂) ≥ ... ≥ Im(λₘ)
  - Allocated by the subroutine

- `eigenvectors` (COMPLEX(rk), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT))
  - Ritz eigenvectors (size n×m)
  - Column i corresponds to eigenvalues(i)
  - Allocated by the subroutine

- `ERROR_STATUS` (INTEGER(ik), INTENT(OUT))
  - `0`: Success
  - `801`: Initial vector has zero norm
  - `802`: Invalid system dimension
  - `803`: Invalid Krylov size
  - `804`: Memory allocation failed
  - `805`: LAPACK eigenvalue computation failed

## Algorithm

1. **Arnoldi Iteration**:
   - Builds orthonormal Krylov basis V = [v₁, v₂, ..., vₘ₊₁]
   - Constructs upper Hessenberg matrix H ((m+1)×m)
   - Uses Gram-Schmidt orthogonalization

2. **Eigenvalue Computation**:
   - Extracts m×m upper block H_m from H
   - Calls LAPACK `DGEEV` to compute eigenvalues and eigenvectors of H_m
   - These are the Ritz values (approximate eigenvalues of A)

3. **Ritz Vector Computation**:
   - Transforms eigenvectors back to full space: Ritz vectors = V * y
   - V is n×m Krylov basis, y are m×m eigenvectors of H_m

4. **Sorting**:
   - Sorts eigenvalues and eigenvectors by descending Im(λ)
   - Most unstable modes (largest positive growth rates) appear first

## Test Matrix

Currently uses a symmetric tridiagonal test matrix:
```
A(i,i) = 0
A(i,i+1) = 1  (superdiagonal)
A(i,i-1) = 1  (subdiagonal)
```

This matrix has known eigenvalues: λₖ = 2cos(kπ/(n+1)) for k=1,...,n

**Note**: This will be replaced with CFD solver calls in the future.

## Configuration

Add to `inputs/inputs.in`:

```fortran
&Arnoldi
    krylov_size = 10,
    frechet_order = 'First',
    eps_0 = 1.0E-6,
    TTime = 1.0,
/
```

## Dependencies

**Modules:**
- `accuracy`: Precision types (ik, rk)
- `error_handling`: Error codes and logging
- `variables`: Global configuration variables

**External Libraries:**
- **LAPACK**: `DGEEV` routine for eigenvalue computation
  - Link with `-llapack -lblas` or use optimized BLAS (MKL, OpenBLAS, Accelerate)

## Compilation

The module requires LAPACK. Update the Makefile to link LAPACK:

```makefile
# For macOS with Accelerate framework:
LIBBLAS := -framework Accelerate

# For Linux with system LAPACK:
# LIBBLAS := -llapack -lblas

# Add to linking step:
$(BINDIR)/$(EXECUT): $(OBJPROG)
	$(FC) $^ $(LIBBLAS) -o $@
```

## Testing

### Standalone Test

Compile and run the test program:

```bash
# From bin/ directory
gfortran -c ../src/accuracy.f90 -o ../obj/accuracy.o -J ../mod/
gfortran -c ../src/error_handling.f90 -o ../obj/error_handling.o -J ../mod/
gfortran -c ../src/variables.f90 -o ../obj/variables.o -J ../mod/
gfortran -c ../src/Arnoldi.f90 -o ../obj/Arnoldi.o -J ../mod/
gfortran -c ../src/test_arnoldi.f90 -o ../obj/test_arnoldi.o -J ../mod/

# Link with LAPACK (macOS):
gfortran ../obj/accuracy.o ../obj/error_handling.o ../obj/variables.o \
         ../obj/Arnoldi.o ../obj/test_arnoldi.o \
         -framework Accelerate -o test_arnoldi

# Run:
./test_arnoldi
```

### Expected Output

```
===============================================
ARNOLDI MODULE - EIGENVALUE COMPUTATION TEST
===============================================

System dimension n = 10
Krylov size m = 5

Generating random initial disturbance vector...

Running Arnoldi iteration with eigenvalue computation...

Arnoldi: Using test matrix A (10x10) - tridiagonal structure

===============================================
SUCCESS: Eigenvalue computation completed!
===============================================

Ritz eigenvalues (sorted by descending Im part):
-----------------------------------------------
  #    Real Part         Imag Part         |λ|
-----------------------------------------------
  1   1.234567E+00      5.678901E-01      1.356789E+00
  ...
```

## Usage Example

```fortran
PROGRAM example
    USE accuracy
    USE Arnoldi
    
    REAL(rk), DIMENSION(100) :: v_init
    COMPLEX(rk), DIMENSION(:), ALLOCATABLE :: evals
    COMPLEX(rk), DIMENSION(:,:), ALLOCATABLE :: evecs
    INTEGER(ik) :: ierr
    
    ! Initialize disturbance vector
    CALL RANDOM_NUMBER(v_init)
    
    ! Compute eigenvalues
    CALL arnoldi_eigenvalues(v_init, 10, 'First', 1.0E-6_rk, 1.0_rk, &
                            evals, evecs, ierr)
    
    IF (ierr /= 0) STOP 'Arnoldi failed'
    
    ! Use results...
    PRINT *, 'Leading eigenvalue:', evals(1)
    
    ! Cleanup
    DEALLOCATE(evals, evecs)
END PROGRAM
```

## Future Integration

The matrix-vector product (currently `MATMUL(A, v)`) will be replaced with:

1. **Write perturbation**: Apply v to flow field
2. **Run CFD solver**: Advance time by ΔT or compute linearized operator
3. **Read result**: Extract resulting perturbation
4. **Return**: w = result of operator action

This allows computing eigenvalues of the linearized Navier-Stokes operator without forming the full matrix.

## Error Handling

All errors are logged through the centralized error handling system:

```fortran
IF (ERROR_STATUS /= 0) THEN
    CALL log_error(ERROR_STATUS)
    ! Handle error appropriately
END IF
```

Error messages include module name, error code, and detailed description.

## Validation

The test matrix eigenvalues can be verified analytically:
- For n=10 tridiagonal matrix: λₖ = 2cos(kπ/11)
- λ₁ ≈ 1.9563, λ₂ ≈ 1.7321, λ₃ ≈ 1.3737, etc.

Compare computed Ritz values against these for validation.

## Status

✅ **Complete**: Full eigenvalue computation with LAPACK
✅ **Complete**: Ritz vector computation
✅ **Complete**: Sorting by descending Im(λ)
✅ **Complete**: Configuration system integration
✅ **Complete**: Error handling
✅ **Complete**: Test program with hardcoded matrix
🔜 **Next**: Replace matrix-vector product with CFD solver calls
🔜 **Next**: Integrate with main TStep workflow
