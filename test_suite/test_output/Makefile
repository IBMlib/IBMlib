#############################################################3
#         Makefile for task.a
# 
############################################################

task.a: test_output.o
	ar rcs task.a test_output.o

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

input_parser.mod  time_tools.mod  run_context.mod  physical_fields.mod  particles.mod output.mod netcdf.mod netcdf_writer.mod ascii_writer.mod: 
	@echo "update of" $@ "suppressed by test_output Makefile"

