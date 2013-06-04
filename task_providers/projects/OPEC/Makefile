########################################################
#                Makefile for task.a                   #
#                                                      # 
#      Set gmake variable MAIN to the desired          #
#      task_provider test (in file $(MAIN).f)          #
########################################################

MAIN = opec1

task.a: $(MAIN).o
	ar rcs task.a $(MAIN).o

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

input_parser.mod  time_tools.mod  run_context.mod  physical_fields.mod: 
	echo "update of" $@ "suppressed"

