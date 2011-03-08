############################################################3
#              Makefile for physical_fields.mod 
#     $Rev: $
#     $LastChangedDate:  $
#     $LastChangedBy: $ 
#
#     Configurations: 
#       Rev 286-310 are based on grid data engine
#       ../generic_elements/regular_lonlat_grid.f
#       Rev 311- are based on the mesh_grid generalization


############################################################
LOCAL_OBJ    = horizontal_grid_transformations.o horizontal_representation.o mesh_grid.o time_services.o MyOcean_daily_mean.o
GENERIC_MESH = ../../generic_elements

physical_fields.a: dependences.mk $(LOCAL_OBJ) 
	ar rcs physical_fields.a $(LOCAL_OBJ)


physical_fields.mod:


# --- interface assembling rules: point to variants that provides sub interfaces

horizontal_grid_transformations.f:
	ln -sf $(GENERIC_MESH)/horizontal_grid_transformations/horizontal_grid_transformations_lonlatgrid.f horizontal_grid_transformations.f

horizontal_representation.f:
	ln -sf  $(GENERIC_MESH)/horizontal_representation/horizontal_representation_scan_staywet.f horizontal_representation.f

time_services.f:
	ln -sf $(GENERIC_MESH)/time_services.f

mesh_grid.f:
	ln -sf $(GENERIC_MESH)/mesh_grid.f  # same name



include $(COMMON_RULES) # shared information between makefiles
FCFLAGS += -I/usr/local/include   # should contain netcdf.mod
VPATH   += /usr/local/include     # let gmake pick up prerequisite netcdf.mod

clean: FORCE
	rm -f  *.o  *.a  *.mod  dependences.mk

FORCE:

# --- Locally tag external modules that should not be build starting from this
#     makefile, but included via the module-search-path-include 
#     compiler directive (ifort: -I). Search paths inherited from 
#     calling makefiles
#     EXTRES = -x time_tools -x run_context

time_tools.mod geometry.mod array_tools.mod constants.mod run_context.mod  input_parser.mod  netcdf.mod:
	@echo "update of" $@ "suppressed"
