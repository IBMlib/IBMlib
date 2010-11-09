##################################################################
#              M a i n   c o n f i g u r a t i o n               #
##################################################################
EXECUTABLE          = tracker   # name of the executable to build
PHYSICAL_FIELDS_DIR = $(IBMLIB_DIR)/oceanography_providers/vortex1
PARTICLE_STATE_DIR  = $(IBMLIB_DIR)/biology_providers/passive
TASK_DIR            = $(IBMLIB_DIR)/task_providers/basic_simulation
##################################################################

#Compiler settings - use defaults as basis
include $(BUILD_TOOLS)/compiler_defaults.mk

#Modify compiler settings 
FCFLAGS += -check bounds -check uninit -traceback  #Add bounds checking and tracebacks
