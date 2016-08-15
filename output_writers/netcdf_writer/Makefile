############################################################3
#                Makefile for netcdf_writer module
############################################################

netcdf_writer.a: netcdf_writer.o 
	ar rcs netcdf_writer.a netcdf_writer.o

netcdf_writer.mod netcdf_writer.o: netcdf_writer.f particle_state.a
	$(FC) -c $(FCFLAGS) $(FPPFLAGS) $< -o $*.o


include $(COMMON_RULES)
#FCFLAGS += -I/usr/local/include   # should contain netcdf.mod

clean: FORCE
	rm -f  *.o  *.a  *.mod  dependences.mk

FORCE:

#Setup generic pseudonyms
archive: netcdf_writer.a 

module: netcdf_writer.mod
       

# --- Locally tag external modules that should not be build from this
#     makefile, but included via the module-search-path-include 
#     compiler directive (ifort: -I). Search paths inherited from 
#     calling makefiles
#
#     EXTRES = -x time_tools  -x particle_tracking  -x physical_fields  -x run_context

particles.mod output.mod netcdf.mod:
	@echo "update of" $@ "suppressed by netcdf_writer makefile"

