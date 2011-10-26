############################################################################
#           Makefile for particle_state.mod for sprat
#     
#     This makefile compiles either multi stage or only feeding larvae setup
#     This is selectedby setting LOCAL_OBJ below
#
#     $Rev:  $
#     $LastChangedDate: $
#     $LastChangedBy: $ 
############################################################################
ONLY_LARVAE = feeding_larval_stage_2_particle_state.o  numerical_1d_integrals.o  particle_state_base.o  prey_community.o  larval_properties.o  feeding_larval_stage.o 

MULTI_STAGE = numerical_1d_integrals.o  particle_state_base.o  prey_community.o  larval_properties.o  feeding_larval_stage.o  yolksac_properties.o  egg_properties.o  released_egg_stage.o  yolksac_larval_stage.o  early_life_stages.o

#
# >>>>>>>> select which variant you want to build by setting LOCAL_OBJ <<<<<<<<<<<
# 
LOCAL_OBJ   = $(ONLY_LARVAE)  

SOURCES     = $(addsuffix .f, $(basename $(LOCAL_OBJ)))

GENERIC_BIOENERGETICS = ../generic_bioenergetics

particle_state.a: $(SOURCES) dependences.mk $(LOCAL_OBJ) 
	ar rcs particle_state.a $(LOCAL_OBJ) 

particle_state.mod:

# ---------------------------------------------------------------------------------
# --- interface assembling rules: point to variants that provides sub interfaces
# ---------------------------------------------------------------------------------

# --- specific components ---

larval_properties.f:
	ln -sf sprat_feeding_larv.f larval_properties.f
yolksac_properties.f:
	ln -sf sprat_yolksac_larv.f yolksac_properties.f
egg_properties.f:
	ln -sf sprat_egg.f          egg_properties.f



# --- generic components ---


early_life_stages.f:
	ln -sf $(GENERIC_BIOENERGETICS)/early_life_stages.f
released_egg_stage.f:
	ln -sf $(GENERIC_BIOENERGETICS)/released_egg_stage.f
yolksac_larval_stage.f:
	ln -sf $(GENERIC_BIOENERGETICS)/yolksac_larval_stage.f
feeding_larval_stage_2_particle_state.f:
	ln -sf $(GENERIC_BIOENERGETICS)/feeding_larval_stage_2_particle_state.f
particle_state_base.f:
	ln -sf $(GENERIC_BIOENERGETICS)/particle_state_base.f
numerical_1d_integrals.f:
	ln -sf $(GENERIC_BIOENERGETICS)/numerical_1d_integrals.f
prey_community.f:
	ln -sf $(GENERIC_BIOENERGETICS)/prey_community.f
feeding_larval_stage.f:
	ln -sf  $(GENERIC_BIOENERGETICS)/optimal_forager.f  feeding_larval_stage.f


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
