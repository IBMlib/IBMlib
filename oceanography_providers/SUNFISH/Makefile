############################################################3
#              Makefile for physical_fields.mod 
############################################################


physical_fields.a: dependences.mk regular_lonlat_grid.o sunfish_pbi.o 
	ar rcs physical_fields.a sunfish_pbi.o regular_lonlat_grid.o
	cd $(IBMLIB_DIR); ln -sf $(shell pwd)/regular_lonlat_grid.mod # tmp fix: project up regular_lonlat_grid.mod

physical_fields.mod:

regular_lonlat_grid.f:
	ln -sf ../generic_elements/regular_lonlat_grid.f


include $(COMMON_RULES) # shared information between makefiles
FCFLAGS += -I/usr/local/include   # should contain netcdf.mod
VPATH   += /usr/local/include     # let gmake pick up prerequisite netcdf.mod

clean: FORCE
	rm -f  *.o  *.a  *.mod  dependences.mk

FORCE:

# --- Locally tag external modules that should not be build starting from this
#     makefile, but included via the module-search-path-include 
#     compiler directive (ifort: -I). Search paths inherited from 
#     calling makefiles
#     EXTRES = -x time_tools -x run_context

time_tools.mod constants.mod run_context.mod  input_parser.mod  netcdf.mod:
	@echo "update of" $@ "suppressed"


