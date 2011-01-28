##################################################################
#              M a i n   c o n f i g u r a t i o n               #
##################################################################
EXECUTABLE          = ibmrun   # name of the executable to build
PHYSICAL_FIELDS_DIR = $(IBMLIB_DIR)/oceanography_providers/linear_field
PARTICLE_STATE_DIR  = $(IBMLIB_DIR)/biology_providers/generic_larvae
TASK_DIR            = $(IBMLIB_DIR)/task_providers/trajectories
OUTPUT_WRITER_DIRS = $(IBMLIB_DIR)/output_writers/ascii_writer
OUTPUT_WRITER_DIRS  += $(IBMLIB_DIR)/output_writers/netcdf_writer
##################################################################

