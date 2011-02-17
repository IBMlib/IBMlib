# -----------------------------------------------------------------
# config.mk should define these 4 targets:
#
#    EXECUTABLE          : name of the final executable
#    PHYSICAL_FIELDS_DIR : directory of the PHYSICAL_FIELDS interface
#    PARTICLE_STATE_DIR  : directory of the PARTICLE_STAT interface
#    TASK_DIR            : directory of the TASK interface
#
#    (each directory above referenced to the base library)
# -----------------------------------------------------------------
EXECUTABLE          = ibmrun
#
PHYSICAL_FIELDS_DIR = $(IBMLIB_DIR)/oceanography_providers/linear_field
#PHYSICAL_FIELDS_DIR = oceanography_providers/POLCOMS+ERSEM
#PHYSICAL_FIELDS_DIR = oceanography_providers/SUNFISH
#
PARTICLE_STATE_DIR  = $(IBMLIB_DIR)/biology_providers/passive
#PARTICLE_STATE_DIR  = biology_providers/sandeel_simple
# 
#TASK_DIR            = task_providers/test
TASK_DIR            = $(IBMLIB_DIR)/task_providers/basic_simulation
#
