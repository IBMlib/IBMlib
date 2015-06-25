############################################################3
#                Makefile for particle_state.mod
############################################################
# 
LOCAL_OBJ = sea_grass_0.o  settler.o 

SOURCES     = $(addsuffix .f, $(basename $(LOCAL_OBJ)))

particle_state.a: $(SOURCES) dependences.mk $(LOCAL_OBJ) 
	ar rcs particle_state.a $(LOCAL_OBJ) 

particle_state.mod:
	

include $(COMMON_RULES)

settler.f:
	ln -sf ../settler_simple/settler.f

clean: FORCE
	rm -f  *.o  *.a  *.mod  dependences.mk

FORCE:

# --- Locally tag external modules that should not be build from this
#     makefile, but included via the module-search-path-include 
#     compiler directive (ifort: -I). Search paths inherited from 
#     calling makefiles
#
#     EXTRES = -x time_tools  -x particle_tracking  -x physical_fields  -x run_context

time_tools.mod  polygons.mod  particle_tracking.mod  physical_fields.mod  run_context.mod  output.mod:
	echo "update of" $@ "suppressed"
