############################################################
#                Makefile for particle_state.mod
############################################################

particle_state.a: particle_state.o
	ar rcs particle_state.a particle_state.o

particle_state.mod particle_state.o: letcher_sandeel.f
	$(FC) -c $(FCFLAGS) $(FPPFLAGS) $< -o $*.o



include $(COMMON_RULES)  # shared information between makefiles


clean: FORCE
	rm -f  *.o  *.a  *.mod  dependences.mk

FORCE:

# --- Locally tag external modules that should not be build from this
#     makefile, but included via the module-search-path-include 
#     compiler directive (ifort: -I). Search paths inherited from 
#     calling makefiles
#
#     EXTRES = -x time_tools  -x particle_tracking  -x physical_fields  -x run_context

time_tools.mod  particle_tracking.mod  physical_fields.mod  run_context.mod:
	echo "update of" $@ "suppressed"
