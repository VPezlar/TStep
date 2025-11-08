# Requires GNU make 3.82 or older. Indentation will not work correctly with other distributions.

SRCDIR := ../src
MODDIR := ../mod
OBJDIR := ../obj
BINDIR := ../bin
EXECUT := TStep

# LIBBLAS := -framework Accelerate
# LIBHDF5 := -I$(shell brew --prefix hdf5)/include -L$(shell brew --prefix hdf5)/lib -lhdf5 -lhdf5_fortran

FC := gfortran
FFLAGS := -O3 -J$(MODDIR) -c -fcheck=all -Warray-bounds -ffree-form -Wall -Wno-tabs #-fno-backtrace

OBJECTS = accuracy.o \
	  error_handling.o \
	  variables.o \
	  setup.o \
	  Arnoldi.o \
	  OpenFOAM_IO.o \
	  read_flow.o \
	  write_output.o \
	  call_CFD.o \
	  random_disturbance.o \
	  write_flow.o \
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
