ROOTDIR = ../../..
BDIR = build
PLIB = pepc-vlib

include $(ROOTDIR)/makefile.defs

all: $(PLIB)

$(PLIB):
	rm -rf $(BDIR); mkdir -p $(BDIR)
	rm -rf $(PLIB).a
	cd $(ROOTDIR)/; $(MAKE) -i $(PLIB) 
	cp $(ROOTDIR)/build/$(MACH)/$(PLIB)/*.f90 $(BDIR)/. 
	cp $(ROOTDIR)/build/$(MACH)/$(PLIB)/*.c $(BDIR)/.
	cp $(ROOTDIR)/build/$(MACH)/$(PLIB)/*.h $(BDIR)/.
	cp $(ROOTDIR)/build/$(MACH)/$(PLIB)/*.a $(BDIR)/.	
	cp $(ROOTDIR)/build/$(MACH)/$(PLIB)/makefile.* $(BDIR)/.
	cp makefile.template $(BDIR)/makefile.build
	cd $(BDIR); $(MAKE) -f makefile.build TARGET=$(PLIB).a
	
clean:
	rm -rf *.o *.mod $(BDIR)

# May remove the rest if you just want the libfile, the following is just for testing 

f90sources   = physvars.f90 files.f90 manipulate_particles.f90 diagnostics.f90 pepc.f90
OBJ = $(f90sources:.f90=.o)
	
pepc-test: $(PLIB) $(OBJ)
	$(LD) $(LDFLAGS) $(OBJ) -o $@ $(LDFLAGS) -I$(BDIR) $(BDIR)/$(PLIB).a $(BDIR)/libsl.a
	
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -I$(BDIR) $(BDIR)/$(PLIB).a $(BDIR)/libsl.a

	

