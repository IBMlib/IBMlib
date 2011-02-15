############################################################
#          Makefile for physical_fields for SUBFISH PBI
#     $Rev$
#     $LastChangedDate$
#     $LastChangedBy$ 
#
#     Configurations: 
#       Rev 273-284 are based on grid data engine
#       ../generic_elements/regular_lonlat_grid.f
#       Rev 285- are based on the mesh_grid generalization
#
############################################################
LOCAL_OBJ = horizontal_grid_transformations.o horizontal_representation.o mesh_grid.o sunfish_pbi.o 

physical_fields.a: dependences.mk $(LOCAL_OBJ) 
	ar rcs physical_fields.a $(LOCAL_OBJ)


physical_fields.mod:


# --- interface assembling rules: point to variants that provides sub interfaces

horizontal_grid_transformations.f:
	ln -sf ../generic_elements/horizontal_grid_transformations/horizontal_grid_transformations_lonlatgrid.f horizontal_grid_transformations.f

horizontal_representation.f:
	ln -sf  ../generic_elements/horizontal_representation/horizontal_representation_staywet.f horizontal_representation.f

mesh_grid.f:
	ln -sf ../generic_elements/mesh_grid.f  # same name



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


