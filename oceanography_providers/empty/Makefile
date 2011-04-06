############################################################3
#              Makefile for physical_fields.mod 
############################################################



physical_fields.a: empty.o
	ar rcs physical_fields.a empty.o

physical_fields.mod physical_fields.o: empty.f
	$(FC) -c $(FCFLAGS) $(FPPFLAGS) $< -o $*.o


include $(COMMON_RULES)

clean: FORCE
	rm -f  *.o  *.a  *.mod  dependences.mk

FORCE:

# --- Locally tag external modules that should not be build from this
#     makefile, but included via the module-search-path-include 
#     compiler directive (ifort: -I). Search paths inherited from 
#     calling makefiles
#     EXTRES = -x time_tools -x run_context

time_tools.mod : 
	echo "update of" $@ "suppressed"


