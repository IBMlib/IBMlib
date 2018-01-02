########################################################
#                Makefile for task.a                   #
#                                                      # 
#      Set gmake variable MAIN to the desired          #
#      task_provider test (in file $(MAIN).f)          #
########################################################

MAIN = connectivity_simulation_ULF

task.a: connectivity.o $(MAIN).o 
	ar rcs task.a connectivity.o  $(MAIN).o

include $(COMMON_RULES)


connectivity.f:
	ln -sf ../generic_elements/connectivity.f

clean: FORCE
	rm -f  *.o  *.a  *.mod  dependences.mk

FORCE:

# --- Locally tag external modules that should not be build from this
#     makefile, but included via the module-search-path-include 
#     compiler directive (ifort: -I). Search paths inherited from 
#     calling makefiles
#
#     EXTRES = -x time_tools  -x run_context  -x physical_fields  -x particles

input_parser.mod  time_tools.mod  run_context.mod  physical_fields.mod  particle_state.mod  netcdf.mod: 
	echo "update of" $@ "suppressed"

