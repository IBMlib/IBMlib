############################################################3
#         Makefile for task.a
# 
############################################################

task.a: part_interface_test.o
	ar rcs task.a part_interface_test.o

include $(COMMON_RULES)

clean: FORCE
	rm -f  *.o  *.a  *.mod  dependences.mk

FORCE:

# --- Locally tag external modules that should not be build from this
#     makefile, but included via the module-search-path-include 
#     compiler directive (ifort: -I). Search paths inherited from 
#     calling makefiles
#
#     EXTRES = -x time_tools  -x run_context  -x physical_fields  -x particles

input_parser.mod  time_tools.mod  run_context.mod  particle_state.mod  particle_tracking.mod  output.mod:
	echo "update of" $@ "suppressed"

