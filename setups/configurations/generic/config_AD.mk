# -----------------------------------------------------------------
#   The config file for an Eulerian configuration is pretty similar
#   to that of a particle setup, apart from PARTICLE_STATE_DIR being 
#   replaced by EULERIAN_DYNAMICS_DIR
#   Please notice that the config file should have the name "config.mk"
#   and be placed in the IBMLIB root directory
#
#    EXECUTABLE            : name of the final executable
#    PHYSICAL_FIELDS_DIR   : directory of the PHYSICAL_FIELDS interface
#    EULERIAN_DYNAMICS_DIR : directory of the EULERIAN interface
#    TASK_DIR              : directory of the TASK interface
#
#    (each directory above referenced to the base library)
# -----------------------------------------------------------------
EXECUTABLE            = adrun
PHYSICAL_FIELDS_DIR   = $(IBMLIB_DIR)/oceanography_providers/SUNFISH_long_phys
EULERIAN_DYNAMICS_DIR = $(IBMLIB_DIR)/eulerian_providers/demo
TASK_DIR              = $(IBMLIB_DIR)/task_providers/eulerian_dynamics
