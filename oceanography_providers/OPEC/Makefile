##########################################################################
#     Makefile for physical_fields for OPEC cmod/ERGOM PBI
#  
#
#     $Rev: 320 $
#     $LastChangedDate: 2011-03-08 11:26:42 +0100 (Tue, 08 Mar 2011) $
#     $LastChangedBy: asch $     
#
#     TODO: add big_endian CF option
##########################################################################
LOCAL_OBJ = horizontal_grid_transformations.o horizontal_representation.o mesh_grid.o time_services.o read_cmod_ergom.o opec.o

GENERIC_MESH = ../generic_elements

physical_fields.a: dependences.mk $(LOCAL_OBJ) 
	ar rcs physical_fields.a $(LOCAL_OBJ)


physical_fields.mod:

# read_cmod.f must be compiled with big_endian option
read_cmod.o: read_cmod.f
	$(FC) -c $(FCFLAGS) $(BIG_ENDIAN) $(FPPFLAGS) $< -o $@


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


