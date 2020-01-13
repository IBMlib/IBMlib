##########################################################################################
#     Makefile for physical_fields for OPEC cmod/ERGOM PBI (two field seamless overlay)
#  
#     New setup based on classes lonlat_grid/lonlat_data2D/lonlat_data3D
#
#
#     $Rev: 320 $
#     $LastChangedDate: 2011-03-08 11:26:42 +0100 (Tue, 08 Mar 2011) $
#     $LastChangedBy: asch $     
#
#     TODO: add biogeochemistry (remove plug, add parallel objects)
##########################################################################################
FCFLAGS     += $(BIG_ENDIAN) # add local compiler flag, defined in compiler_defaults.mk

#
LOCAL_OBJ = time_services.o read_cmod_ergom.o regular_lonlat_grid.o regular_lonlat_data.o opec_two_set_fusion.o

GENERIC_MESH   = ../generic_elements
GENERIC_LONLAT = ../generic_elements/regular_lonlat

physical_fields.a: dependences.mk $(LOCAL_OBJ) 
	ar rcs physical_fields.a $(LOCAL_OBJ)

physical_fields.mod:


# --- interface assembling rules: point to variants that provides sub interfaces

time_services.f:
	ln -sf $(GENERIC_MESH)/time_services.f

read_cmod_ergom.f:
	ln -sf ../OPEC/read_cmod_ergom.f


regular_lonlat_grid.f: coastline_intersection_for_grid.f
	ln -sf $(GENERIC_LONLAT)/regular_lonlat_grid.f

regular_lonlat_data.f:
	ln -sf $(GENERIC_LONLAT)/regular_lonlat_data.f

# ---- link to particular implementation of subroutine coastline_intersection_for_grid() for inline in regular_lonlat_grid.f
coastline_intersection_for_grid.f:
	ln -sf  $(GENERIC_LONLAT)/coastline_intersection_for_grid_staywet.f  coastline_intersection_for_grid.f

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

time_tools.mod geometry.mod array_tools.mod constants.mod run_context.mod  input_parser.mod:
	@echo "update of" $@ "suppressed"


