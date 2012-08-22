#############################################################################
#             Makefile for advection_diffusion.mod
# 
#  Assemble and compile a provider of the advection_diffusion interface
#############################################################################

advection_diffusion.a: advection_diffusion.o
	ar rcs advection_diffusion.a  eulerian_production_rate.o advection_diffusion.o

advection_diffusion.mod advection_diffusion.o: advection_diffusion.f eulerian_production_rate.f
	$(FC) -c $(FCFLAGS) $(FPPFLAGS) $< -o $*.o


advection_diffusion.f:
	ln -sf ../generic_advection_diffusion/advection_diffusion_mesh_grid.f         advection_diffusion.f

eulerian_production_rate.f:
	ln -sf ../generic_eulerian_production_rate/logistic_with_spatial_mortality.f  eulerian_production_rate.f


include $(COMMON_RULES)


clean: FORCE
	rm -f  *.o  *.a  *.mod  dependences.mk

FORCE:

# --- Locally tag external modules that should not be build from this
#     makefile, but included via the module-search-path-include 
#     compiler directive (ifort: -I). Search paths inherited from 
#     calling makefiles 
#
#     EXTRES = -x time_tools  -x particle_tracking  -x physical_fields  -x run_context

physical_fields.mod  mesh_grid.mod  horizontal_grid_transformations.mod  geometry.mod  run_context.mod  input_parser:
	echo "update of" $@ "suppressed"
