# Arnoldi Module Rewrite - Summary

## Changes Overview

The Arnoldi module has been completely rewritten to compute eigenvalues and eigenvectors with the specifications you provided.

## Files Modified

### 1. **src/Arnoldi.f90** - Complete rewrite
- New subroutine: `arnoldi_eigenvalues()` (replaces `arnoldi_iter`)
- **Inputs**: v_init, m, frechet_order, eps_0, TTime
- **Outputs**: eigenvalues (complex), eigenvectors (complex), ERROR_STATUS
- Uses LAPACK `DGEEV` for eigenvalue computation
- Computes Ritz values (eigenvalues of H_m, the m×m upper block of Hessenberg matrix)
- Computes Ritz vectors (eigenvectors transformed to full space: V*y)
- Sorts by **descending imaginary part**: Im(λ)
- Hardcoded test matrix: tridiagonal with known eigenvalues

### 2. **src/variables.f90** - Added Arnoldi parameters
```fortran
INTEGER(ik) :: krylov_size
CHARACTER(len=20) :: frechet_order
REAL(rk) :: eps_0, TTime
```

### 3. **src/setup.f90** - Added Arnoldi namelist
- New namelist: `&Arnoldi`
- Reads: `krylov_size`, `frechet_order`, `eps_0`, `TTime`
- Integrated into configuration reading workflow

### 4. **src/error_handling.f90** - Added error codes
- `ERR_ARNOLDI_ALLOC = 804`: Memory allocation failed
- `ERR_ARNOLDI_LAPACK = 805`: LAPACK computation failed

### 5. **src/test_arnoldi.f90** - Updated test program
- Tests new eigenvalue API
- Displays sorted eigenvalues and eigenvectors
- Validates against hardcoded test matrix

### 6. **bin/makefile** - Already updated
- Arnoldi.o already in build order
- **NOTE**: Need to uncomment LAPACK linking for compilation

## Key Features

### Eigenvalue Computation
- ✅ Arnoldi iteration builds Krylov basis V and Hessenberg matrix H
- ✅ Extracts m×m upper block H_m
- ✅ LAPACK DGEEV computes eigenvalues and eigenvectors of H_m
- ✅ Returns Ritz values (approximate eigenvalues of A)
- ✅ Returns Ritz vectors (approximate eigenvectors of A)

### Sorting
- ✅ Eigenvalues sorted by **descending Im(λ)**
- ✅ Eigenvectors reordered to match
- ✅ Most unstable modes appear first

### Test Matrix
- ✅ Tridiagonal matrix: A(i,i)=0, A(i,i±1)=1
- ✅ Known eigenvalues: λₖ = 2cos(kπ/(n+1))
- ✅ Size inferred from v_init dimension
- ⚠️ Will be replaced with CFD solver calls

### Error Handling
- ✅ Validates all inputs
- ✅ Checks allocation status
- ✅ Handles LAPACK errors
- ✅ Integrated with centralized error system

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

## Compilation Requirements

**LAPACK is required**. Update `bin/makefile`:

```makefile
# Uncomment one of these:

# macOS:
LIBBLAS := -framework Accelerate

# Linux with system LAPACK:
# LIBBLAS := -llapack -lblas

# Intel MKL:
# LIBBLAS := -lmkl_rt

# Already in makefile (line 36):
$(BINDIR)/$(EXECUT): $(OBJPROG)
	$(FC) $^ $(LIBHDF5) $(LIBBLAS) -o $@
```

## Testing on HPC

### Compile Test Program

```bash
cd bin

# Compile dependencies
gfortran -c ../src/accuracy.f90 -o ../obj/accuracy.o -J ../mod/
gfortran -c ../src/error_handling.f90 -o ../obj/error_handling.o -J ../mod/
gfortran -c ../src/variables.f90 -o ../obj/variables.o -J ../mod/

# Compile Arnoldi module
gfortran -c ../src/Arnoldi.f90 -o ../obj/Arnoldi.o -J ../mod/

# Compile test program
gfortran -c ../src/test_arnoldi.f90 -o ../obj/test_arnoldi.o -J ../mod/

# Link with LAPACK
gfortran ../obj/accuracy.o ../obj/error_handling.o ../obj/variables.o \
         ../obj/Arnoldi.o ../obj/test_arnoldi.o \
         -llapack -lblas -o test_arnoldi

# Or use module system if available:
# module load lapack
# gfortran ... -llapack -lblas -o test_arnoldi

# Run
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
  1   x.xxxxxxE+xx      x.xxxxxxE+xx      x.xxxxxxE+xx
  2   x.xxxxxxE+xx      x.xxxxxxE+xx      x.xxxxxxE+xx
  ...
```

## Validation

The tridiagonal test matrix has known eigenvalues:
- λₖ = 2cos(kπ/(n+1)) for k=1,...,n

For n=10:
- λ₁ ≈ 1.9563
- λ₂ ≈ 1.7321
- λ₃ ≈ 1.3737
- λ₄ ≈ 0.8794
- λ₅ ≈ 0.3473

**Note**: These are real eigenvalues, so Im(λ) ≈ 0 for all. The test matrix eigenvalues are **real**, so sorting by Im(λ) may show zeros or numerical noise.

## Next Steps

1. **Test on HPC**: Compile and run to verify LAPACK integration
2. **Validate Results**: Check computed Ritz values against analytical values
3. **Replace Test Matrix**: Integrate CFD solver calls for matrix-vector products
4. **Main Integration**: Call from main TStep workflow
5. **Output**: Write eigenvalues/eigenvectors to file

## API Summary

```fortran
! New API
CALL arnoldi_eigenvalues(v_init, m, frechet_order, eps_0, TTime, &
                        eigenvalues, eigenvectors, ERROR_STATUS)

! Where:
!   v_init        - REAL(rk), DIMENSION(n)          [IN]
!   m             - INTEGER(ik)                     [IN]
!   frechet_order - CHARACTER(len=*)                [IN]
!   eps_0         - REAL(rk)                        [IN]
!   TTime         - REAL(rk)                        [IN]
!   eigenvalues   - COMPLEX(rk), DIMENSION(m)       [OUT, ALLOCATABLE]
!   eigenvectors  - COMPLEX(rk), DIMENSION(n,m)     [OUT, ALLOCATABLE]
!   ERROR_STATUS  - INTEGER(ik)                     [OUT]
```

## Documentation

See `ARNOLDI_EIGENVALUE_API.md` for complete API documentation.

## Questions?

If you encounter any issues during compilation or testing on HPC:
1. Check LAPACK availability: `module avail lapack` or `ldconfig -p | grep lapack`
2. Verify gfortran version supports BLOCK constructs (>= 4.6)
3. Check error messages for missing symbols (link order matters)
4. Ensure mod files are accessible via `-J` flag
