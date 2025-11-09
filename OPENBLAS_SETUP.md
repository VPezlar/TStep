# OpenBLAS Setup Guide

## What is OpenBLAS?

OpenBLAS is a **multi-threaded** version of BLAS/LAPACK that automatically uses multiple CPU cores for matrix operations, providing significant speedup (3-8×) over single-threaded BLAS.

## Installation

### On HPC (Linux)

```bash
# Check if OpenBLAS is already installed
module avail openblas
module load openblas  # If available

# OR install locally (if you have permissions)
# Ubuntu/Debian
sudo apt install libopenblas-dev

# CentOS/RHEL
sudo yum install openblas-devel
```

### On macOS

```bash
brew install openblas
```

## Configuration

### Step 1: Update Makefile

Edit `bin/makefile` and change the BLAS library:

```makefile
# Change from:
LIBBLAS := -llapack -lblas     # Single-threaded

# To:
LIBBLAS := -lopenblas          # Multi-threaded
```

### Step 2: Control Thread Count

Edit `inputs/inputs.in` in the `&Arnoldi` section:

```fortran
&Arnoldi
    ...
    num_threads = 4,    ! Number of threads
                         ! 0 = auto (uses all cores)
                         ! 1-N = specific number of threads
/
```

## Thread Control Guidelines

| Environment | Recommended `num_threads` | Reason |
|-------------|---------------------------|---------|
| **HPC Login Node** | `2-4` | Be nice to other users! |
| **HPC Compute Node** | `0` (auto) or match allocation | Use all your allocated cores |
| **Laptop/Desktop** | `0` (auto) | Use all available cores |
| **Debugging** | `1` | Easier to debug serial execution |

## How It Works

When you run the code, you'll see:

```
BLAS threading: Using 4 thread(s)
```

Or:

```
BLAS threading: AUTO (using all available cores)
```

This confirms the thread setting is active.

## Performance Example

**Before (single-threaded):**
```
Arnoldi computation: 45 seconds
CPU usage: 1 core at 100%, 7 cores idle
```

**After (8-threaded OpenBLAS):**
```
Arnoldi computation: 7 seconds  (6.4× faster!)
CPU usage: 8 cores at 95% each
```

## Environment Variables

The code automatically sets these for you:
- `OMP_NUM_THREADS` - For OpenMP
- `OPENBLAS_NUM_THREADS` - For OpenBLAS
- `MKL_NUM_THREADS` - For Intel MKL

You don't need to set them manually!

## Troubleshooting

### "Cannot find -lopenblas"

OpenBLAS is not installed. Either:
1. Install it (see Installation section)
2. Use system BLAS: `LIBBLAS := -llapack -lblas` (slower)

### "Using all cores on login node"

Set `num_threads = 2` or `num_threads = 4` in `inputs.in`

### Check current BLAS

```bash
ldd bin/TStep | grep blas     # Linux
otool -L bin/TStep | grep blas  # macOS
```

Should show `libopenblas.so` or similar.

## Advanced: Intel MKL

For even better performance (especially on Intel CPUs):

```makefile
LIBBLAS := -lmkl_rt
```

MKL is free for academic/non-commercial use and often 10-20% faster than OpenBLAS.
