############################################################3
#                Makefile for ascii_writer module
############################################################

ascii_writer.a: ascii_writer.o
	ar rcs ascii_writer.a ascii_writer.o

ascii_writer.mod ascii_writer.o: ascii_writer.f
	$(FC) -c $(FCFLAGS) $(FPPFLAGS) $< -o $*.o



include $(COMMON_RULES)


clean: FORCE
	rm -f  *.o  *.a  *.mod  dependences.mk

FORCE:

#Setup generic pseudonyms
archive: ascii_writer.a 

module: ascii_writer.mod

# --- Locally tag external modules that should not be build from this
#     makefile, but included via the module-search-path-include 
#     compiler directive (ifort: -I). Search paths inherited from 
#     calling makefiles
#
#     EXTRES = -x time_tools  -x particle_tracking  -x physical_fields  -x run_context

particles.mod output.mod :
	@echo "update of" $@ "suppressed in ascii_writer makefile"

