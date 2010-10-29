############################################################3
#              Makefile for physical_fields.mod 
############################################################



physical_fields.a: physical_fields.o
	ar rcs physical_fields.a physical_fields.o

physical_fields.mod physical_fields.o: sunfish_pbi.f
	$(FC) -c $(FCFLAGS) $(FPPFLAGS) $< -o $*.o


include $(IBMLIB_DIR)/common.mk   # shared information between makefiles
FCFLAGS += -I/usr/local/include   # should contain netcdf.mod

clean: FORCE
	rm -f  *.o  *.a  *.mod  dependences.mk

FORCE:

# --- Locally tag external modules that should not be build from this
#     makefile, but included via the module-search-path-include 
#     compiler directive (ifort: -I). Search paths inherited from 
#     calling makefiles
#     EXTRES = -x time_tools -x run_context

time_tools.mod constants.mod run_context.mod  input_parser.mod  netcdf.mod:
	echo "update of" $@ "suppressed"


