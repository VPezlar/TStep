# Requires GNU make 3.82 or older. Indentation will not work correctly with other distributions.

SRCDIR := ../src
MODDIR := ../mod
OBJDIR := ../obj
BINDIR := ../bin
EXECUT := TStep

# LAPACK/BLAS is REQUIRED for Arnoldi module eigenvalue computation
# Uncomment one of the following based on your system:
#
# Option 1: OpenBLAS (RECOMMENDED - automatically multi-threaded)
LIBBLAS := -lopenblasp                                 # Linux/macOS with OpenBLAS
#                                                      # Install: sudo apt install libopenblasp-dev (Ubuntu/Debian)
#                                                      #          sudo yum install openblas-devel (CentOS/RHEL)
#                                                      #          brew install openblas (macOS)
#
# Option 2: Intel MKL (Fastest, auto-threaded, free for non-commercial)
# LIBBLAS := -lmkl_rt                                 # Intel MKL
#
# Option 3: System BLAS (Single-threaded, slowest)
# LIBBLAS := -llapack -lblas                          # Linux with system LAPACK (NOT threaded)
#
# Option 4: Apple Accelerate (macOS only, auto-threaded)
# LIBBLAS := -framework Accelerate                    # macOS Accelerate framework

# LIBHDF5 := -I$(shell brew --prefix hdf5)/include -L$(shell brew --prefix hdf5)/lib -lhdf5 -lhdf5_fortran

FC := gfortran
FFLAGS := -O3 -J$(MODDIR) -c -fcheck=all -Warray-bounds -ffree-form -Wall -Wno-tabs #-fno-backtrace

OBJECTS = accuracy.o \
	  error_handling.o \
	  variables.o \
	  setup.o \
	  state_vector.o \
	  frechet_stencil.o \
	  OpenFOAM_IO.o \
	  read_flow.o \
	  write_flow.o \
	  call_CFD.o \
	  Arnoldi.o \
	  write_eigendata.o \
	  write_output.o \
	  random_disturbance.o \
	  main.o

# Compilation rule
$(OBJDIR)/%.o : $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) $< -o $@

OBJPROG = $(addprefix $(OBJDIR)/, $(OBJECTS))

all : $(BINDIR)/$(EXECUT)

$(BINDIR)/$(EXECUT): $(OBJPROG)
	$(FC) $^ $(LIBHDF5) $(LIBBLAS) -o $@

clean:
	rm -f $(OBJDIR)/*.o $(MODDIR)/*.mod $(BINDIR)/$(EXECUT) ./fort.*
