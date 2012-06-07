############################################################################
#           Makefile for particle_state.mod for sprat with MEECE ecotox
#     
#     This makefile supports only multi stage feeding larvae setup
#
#     $Rev:  $
#     $LastChangedDate: $
#     $LastChangedBy: $ 
############################################################################


MULTI_STAGE = feeding_larval_stage.o  released_egg_stage.o  yolksac_larval_stage.o  early_life_stages.o

#
# >>>>>>>> select which variant you want to build by setting LOCAL_OBJ <<<<<<<<<<<
# 
LOCAL_OBJ   = $(MULTI_STAGE)  

SOURCES     = $(addsuffix .f, $(basename $(LOCAL_OBJ)))

GENERIC_BIOENERGETICS = ../generic_bioenergetics

particle_state.a: $(SOURCES) dependences.mk $(LOCAL_OBJ) 
	ar rcs particle_state.a $(LOCAL_OBJ) 

particle_state.mod:


# --- generic components ---

early_life_stages.f:
	ln -sf $(GENERIC_BIOENERGETICS)/early_life_stages.f
released_egg_stage.f:
	ln -sf sprat_egg.f           released_egg_stage.f
yolksac_larval_stage.f:
	ln -sf sprat_yolksac_larv.f  yolksac_larval_stage.f
feeding_larval_stage.f:
	ln -sf  $(GENERIC_BIOENERGETICS)/passive_feeder.f  feeding_larval_stage.f


include $(COMMON_RULES)


clean: FORCE
	rm -f  *.o  *.a  *.mod  dependences.mk $(SOURCES)   # SOURCES are soft links

FORCE:

# --- Locally tag external modules that should not be build from this
#     makefile, but included via the module-search-path-include 
#     compiler directive (ifort: -I). Search paths inherited from 
#     calling makefiles
#
#     EXTRES = -x time_tools  -x particle_tracking  -x physical_fields  -x run_context

time_tools.mod  particle_tracking.mod  physical_fields.mod  run_context.mod spline.mod output.mod:
	echo "update of" $@ "suppressed"
