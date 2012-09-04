DEBUG    = 0


LIBS     = -I$(FFTW_ROOT)/include -lfftw3 -llapack

ifeq ($(DEBUG), 1)
  FFLAGS   =-O0 -C -check noarg_temp_created -g -openmp
else
  FFLAGS   = -O5 -openmp 
endif

EXECUTABLES = spatially_resolved_ccf spatially_resolved_ccf_jackknife spatially_resolved_spherical_fields spatially_resolved_spherical_fields_jackknife spatially_resolved_spherical_spectrum_jackknife
MODULES     = progress.o module_math.o
FC          = mpif90

all: loadmodules $(EXECUTABLES) bfunc.so

clean:
	rm -rf $(EXECUTABLES) *.mod *.o *.so

%: %.f90 makefile $(MODULES)
	@echo "== compiling F90 file :" $<
	$(FCPRE) $(FC) $(FFLAGS) $< -o $@ $(MODULES) $(LIBS)
	
%.o: %.f90 makefile
	@echo "== compiling F90 file :" $<
	$(FCPRE) $(FC) $(FFLAGS) -fPIC $< -c -o $@ $(LIBS)
	
bfunc.so: bfunc_python.f90 module_math.o
	f2py -m bfunc -c module_math.o bfunc_python.f90

loadmodules:
	@echo "Do not forget to 'source modules.sh' before compiling and running"

# this rule prevents the .o-files from being deleted
.SECONDARY: $(MODULES)
