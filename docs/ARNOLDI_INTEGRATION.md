# Arnoldi Module Integration

## Summary

The `Arnoldi.f90` module has been integrated into the TStep codebase. The module is now ready for standalone testing and can be used independently before full integration into the main workflow.

## Changes Made

### 1. Refactored `Arnoldi.f90` to Match Codebase Standards

**Code Quality Improvements:**
- Added `IMPLICIT NONE` declaration
- Added `PRIVATE` and `PUBLIC` declarations for proper encapsulation
- Converted to uppercase Fortran keywords (consistent with codebase style)
- Added comprehensive subroutine documentation header
- Changed `v` intent from `INTENT(OUT)` to `INTENT(INOUT)` (caller provides v(:,1))
- Added `ERROR_STATUS` output parameter for error handling
- Used intrinsic `NORM2` function instead of manual `sqrt(sum(x**2))`
- Replaced `tiny(1.0_rk)` checks with explicit `1.0E-12_rk` threshold
- Improved warning messages for Arnoldi breakdown
- Added proper cleanup on breakdown (zero out remaining entries)

**Error Handling Integration:**
- Added error codes to `error_handling.f90` module:
  - `ERR_ARNOLDI_ZERO_V1 = 801`: Initial vector has zero norm
  - `ERR_ARNOLDI_INVALID_DIM = 802`: Invalid matrix dimension (reserved)
  - `ERR_ARNOLDI_INVALID_KRYLOV = 803`: Invalid Krylov size (reserved)
- Integrated with centralized `log_error()` function
- Added module ID 8 for Arnoldi (error range 800-899)

### 2. Updated Makefile (`bin/makefile`)
- Added `Arnoldi.o` to the `OBJECTS` list after `setup.o`
- Compilation order ensures dependencies are met:
  - `accuracy.o` compiled first (no dependencies)
  - `error_handling.o` compiled before `Arnoldi.o` (Arnoldi uses error handling)

## Module Dependencies

```
Arnoldi
  ├── accuracy (provides: ik, rk precision types)
  └── error_handling (provides: error codes and log_error function)
```

## Compilation

The module will be compiled as part of the standard build process:

```bash
cd bin
make clean
make
```

On the HPC system, this will create:
- Object file: `obj/Arnoldi.o`
- Module file: `mod/arnoldi.mod`

## Standalone Testing

A test program `test_arnoldi.f90` has been created to verify the module works correctly. This program:
1. Creates a simple 5×5 test matrix
2. Initializes a random starting vector
3. Runs the Arnoldi iteration with m=3 Krylov subspace size
4. Outputs the resulting Hessenberg matrix

To compile and run the test (on HPC):
```bash
# Compile dependencies
gfortran -c ../src/accuracy.f90 -o ../obj/accuracy.o -J ../mod/
gfortran -c ../src/error_handling.f90 -o ../obj/error_handling.o -J ../mod/
gfortran -c ../src/variables.f90 -o ../obj/variables.o -J ../mod/
gfortran -c ../src/setup.f90 -o ../obj/setup.o -J ../mod/
gfortran -c ../src/Arnoldi.f90 -o ../obj/Arnoldi.o -J ../mod/

# Compile test program
gfortran -c ../src/test_arnoldi.f90 -o ../obj/test_arnoldi.o -J ../mod/

# Link
gfortran ../obj/accuracy.o ../obj/error_handling.o ../obj/variables.o ../obj/setup.o ../obj/Arnoldi.o ../obj/test_arnoldi.o -o test_arnoldi

# Run
./test_arnoldi
```

## Module API

### `arnoldi_iter(A, n, m, v, H, ERROR_STATUS)`

Implements the Arnoldi iteration algorithm to build an orthonormal basis of the Krylov subspace K_m(A, v1) = span{v1, Av1, A²v1, ..., A^(m-1)v1}.

**Parameters:**
- `A` (REAL(rk), DIMENSION(n,n), INTENT(IN)): Input matrix (n × n)
- `n` (INTEGER(ik), INTENT(IN)): Matrix dimension
- `m` (INTEGER(ik), INTENT(IN)): Krylov subspace size
- `v` (REAL(rk), DIMENSION(n,m+1), INTENT(INOUT)): Basis vectors
  - **INPUT**: `v(:,1)` must contain a non-zero starting vector
  - **OUTPUT**: All columns contain orthonormal Krylov basis vectors
- `H` (REAL(rk), DIMENSION(m+1,m), INTENT(OUT)): Hessenberg matrix
  - Upper Hessenberg form: H(i,j) = 0 for i > j+1
  - Contains projection coefficients from Gram-Schmidt process
- `ERROR_STATUS` (INTEGER(ik), INTENT(OUT)): Error code
  - `0`: Success
  - `801`: Initial vector v(:,1) has zero or near-zero norm

**Algorithm Steps:**
1. Validates and normalizes initial vector `v(:,1)` to unit length
2. For j = 1 to m:
   - **Step 3**: Compute w = A·v_j (matrix-vector product)
   - **Step 4-6**: Orthogonalize w via Gram-Schmidt against v_1, ..., v_j
   - **Step 7-8**: Compute norm ||w|| and store as H(j+1,j)
   - **Step 9-10**: Normalize w and store as v_{j+1}
3. Handles breakdown gracefully (early termination if Krylov space exhausted)

**Error Handling:**
- Returns `ERR_ARNOLDI_ZERO_V1` if initial vector has ||v1|| < 1.0E-12
- Logs error via centralized `log_error()` function
- On breakdown (||w|| ≈ 0), logs warning and zeros remaining entries

**Numerical Stability:**
- Uses intrinsic `NORM2()` for L2-norm calculation
- Zero threshold: `1.0E-12_rk` for double precision
- Classical Gram-Schmidt orthogonalization (may be replaced with MGS if needed)

**Future Modifications:**
- Step 3 (MATMUL) will be replaced with CFD solver calls
- The operator A will represent the linearized Navier-Stokes operator
- Matrix-free implementation: only matrix-vector products needed

## Future Integration

The Arnoldi module is designed to be integrated into the main time-stepping workflow. Step 3 of the algorithm (matrix-vector multiplication) will be replaced with:
- External CFD solver calls
- Flow field perturbations
- Result collection

## Code Quality Standards

The Arnoldi module now follows all TStep coding conventions:

✅ **Style**: Uppercase Fortran keywords, consistent indentation  
✅ **Error Handling**: Integrated with centralized error system (800-899 range)  
✅ **Documentation**: Comprehensive subroutine headers with algorithm details  
✅ **Encapsulation**: Proper `PRIVATE`/`PUBLIC` declarations  
✅ **Precision**: Uses `rk` and `ik` from `accuracy` module  
✅ **Robustness**: Validates inputs, handles edge cases gracefully  
✅ **Intrinsics**: Uses `NORM2`, `DOT_PRODUCT`, `MATMUL` from Fortran standard  

## Error Code Allocation

**Module Range: 800-899 (Arnoldi)**

| Code | Constant | Description |
|------|----------|-------------|
| 801 | `ERR_ARNOLDI_ZERO_V1` | Initial vector has zero or near-zero norm |
| 802 | `ERR_ARNOLDI_INVALID_DIM` | Invalid matrix dimension (reserved) |
| 803 | `ERR_ARNOLDI_INVALID_KRYLOV` | Invalid Krylov subspace size (reserved) |

## Status

✅ **Complete**: Module refactored to match codebase standards  
✅ **Complete**: Error handling integrated (error codes 800-899)  
✅ **Complete**: Comprehensive documentation added  
✅ **Complete**: Test program updated with error checking  
✅ **Complete**: Added to Makefile with correct dependencies  
🔜 **Next**: Test on HPC system  
🔜 **Next**: Replace Step 3 with CFD solver integration
