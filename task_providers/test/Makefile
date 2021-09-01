########################################################
#                Makefile for task.a                   #
#                                                      # 
#      Set gmake variable MAIN to the desired          #
#      task_provider test (in file $(MAIN).f)          #
########################################################

PROVIDERS = test2_pbi.o write_netcdf_lonlat_data.o
#PROVIDERS = test_pbi.o

task.a: $(PROVIDERS)
	ar rcs task.a $(PROVIDERS)

include $(COMMON_RULES)

write_netcdf_lonlat_data.f:
	ln -sf ../generic_elements/write_netcdf_lonlat_data.f  # same name

clean: FORCE
	rm -f  *.o  *.a  *.mod  dependences.mk

FORCE:

# --- Locally tag external modules that should not be build from this
#     makefile, but included via the module-search-path-include 
#     compiler directive (ifort: -I). Search paths inherited from 
#     calling makefiles
#
#     EXTRES = -x time_tools  -x run_context  -x physical_fields  -x particles

input_parser.mod  time_tools.mod  run_context.mod  physical_fields.mod  particles.mod netcdf.mod: 
	echo "update of" $@ "suppressed"

