# Arnoldi Module Integration

## Summary

The `Arnoldi.f90` module has been integrated into the TStep codebase. The module is now ready for standalone testing and can be used independently before full integration into the main workflow.

## Changes Made

### 1. Fixed Precision References in `Arnoldi.f90`
- **Line 36**: Changed `tiny(1.0_dp)` → `tiny(1.0_rk)`
- **Line 64**: Changed `tiny(1.0_dp)` → `tiny(1.0_rk)`

**Reason**: The `accuracy` module defines `rk` (real64) as the standard real precision type, not `dp`. This ensures consistency with the rest of the codebase.

### 2. Updated Makefile (`bin/makefile`)
- Added `Arnoldi.o` to the `OBJECTS` list after `setup.o`
- Compilation order ensures dependencies are met:
  - `accuracy.o` compiled first (no dependencies)
  - `setup.o` compiled before `Arnoldi.o` (Arnoldi depends on both accuracy and setup)

## Module Dependencies

```
Arnoldi
  ├── accuracy (provides: ik, rk)
  └── setup (provides: get_unit and configuration utilities)
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

### `arnoldi_iter(A, n, m, v, H)`

Implements the Arnoldi iteration algorithm to build an orthonormal basis of the Krylov subspace.

**Parameters:**
- `A` (REAL(rk), dimension(n,n), intent(in)): Input matrix
- `n` (INTEGER(ik), intent(in)): Dimension of the matrix
- `m` (INTEGER(ik), intent(in)): Krylov subspace size
- `v` (REAL(rk), dimension(n,m+1), intent(out)): Output basis vectors
  - **Note**: `v(:,1)` should be initialized with a non-zero starting vector before calling
- `H` (REAL(rk), dimension(m+1,m), intent(out)): Output Hessenberg matrix

**Algorithm Steps:**
1. Normalizes the initial vector `v(:,1)` (must be provided)
2. For j = 1 to m:
   - Compute w = A·v_j (Step 3 - will be replaced later)
   - Orthogonalize w against all previous v_i vectors (Gram-Schmidt)
   - Store coefficients in H(i,j)
   - Normalize and store as v_{j+1}
3. Handles breakdown gracefully when ||w|| ≈ 0

**Notes:**
- Step 3 (matrix-vector multiplication) is currently implemented as `w = matmul(A, v(:,j))`
- This step will be replaced in the future with calls to the CFD solver
- The module currently assumes `v(:,1)` is provided and non-zero
- Uses robust zero-checking with `tiny(1.0_rk)` for numerical stability

## Future Integration

The Arnoldi module is designed to be integrated into the main time-stepping workflow. Step 3 of the algorithm (matrix-vector multiplication) will be replaced with:
- External CFD solver calls
- Flow field perturbations
- Result collection

## Status

✅ **Complete**: Module compiles with correct dependencies  
✅ **Complete**: Precision types fixed (`rk` instead of `dp`)  
✅ **Complete**: Added to Makefile  
✅ **Complete**: Test program created  
🔜 **Next**: Test on HPC system  
🔜 **Next**: Replace Step 3 with CFD solver integration
