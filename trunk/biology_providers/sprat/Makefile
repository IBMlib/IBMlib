############################################################3
#                Makefile for particle_state.mod
#
#     $Rev:  $
#     $LastChangedDate: $
#     $LastChangedBy: $ 
############################################################
LOCAL_OBJ = particle_state_base.o prey_community.o sprat.o optimal_forager.o 

GENERIC_BIOENERGETICS = ../generic_bioenergetics

particle_state.a: dependences.mk $(LOCAL_OBJ) 
	ar rcs particle_state.a $(LOCAL_OBJ) 

particle_state.mod:


# --- interface assembling rules: point to variants that provides sub interfaces

particle_state_base.f:
	ln -sf $(GENERIC_BIOENERGETICS)/particle_state_base.f

prey_community.f:
	ln -sf $(GENERIC_BIOENERGETICS)/prey_community.f

optimal_forager.f:
	ln -sf $(GENERIC_BIOENERGETICS)/optimal_forager.f


include $(COMMON_RULES)


clean: FORCE
	rm -f  *.o  *.a  *.mod  dependences.mk

FORCE:

# --- Locally tag external modules that should not be build from this
#     makefile, but included via the module-search-path-include 
#     compiler directive (ifort: -I). Search paths inherited from 
#     calling makefiles
#
#     EXTRES = -x time_tools  -x particle_tracking  -x physical_fields  -x run_context

time_tools.mod  particle_tracking.mod  physical_fields.mod  run_context.mod spline.mod output.mod:
	echo "update of" $@ "suppressed"
