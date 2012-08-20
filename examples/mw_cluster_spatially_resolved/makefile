DEBUG    = 0


LIBS     = -I$(FFTW_ROOT)/include -lfftw3 -llapack

ifeq ($(DEBUG), 1)
  FFLAGS   =-O0 -C -check noarg_temp_created -g -openmp
else
  FFLAGS   = -O5 -openmp 
endif

EXECUTABLES = spatially_resolved_ccf spatially_resolved_ccf_jackknife spatially_resolved_spherical_fields
FC          = mpif90

all: loadmodules $(EXECUTABLES)

clean:
	rm -rf $(EXECUTABLES) *.mod

%: %.f90
	@echo "== compiling F90 file :" $<
	$(FCPRE) $(FC) -c $(FFLAGS) $< -o $@ $(LIBS)

loadmodules:
	@echo "Do not forget to 'source modules.sh' before compiling and running"
	@echo $(FFTW_ROOT)
