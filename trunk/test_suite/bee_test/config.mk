##################################################################
#              M a i n   c o n f i g u r a t i o n               #
##################################################################
EXECUTABLE          = ibmrun_test   # name of the executable to build
PHYSICAL_FIELDS_DIR = $(IBMLIB_DIR)/oceanography_providers/linear_field
PARTICLE_STATE_DIR  = $(IBMLIB_DIR)/biology_providers/passive
TASK_DIR            = $(IBMLIB_DIR)/test_suite/bee_test
##################################################################

#Modify compiler settings 
#bounds checking and tracebacks for ifort are included in compiler_defaults_safe.mk
#FCFLAGS += -check bounds -check uninit -traceback  #Add bounds checking and tracebacks
