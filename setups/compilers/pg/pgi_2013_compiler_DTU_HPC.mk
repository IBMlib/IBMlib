#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
# ---------------------------------------------------
# Portland setup for HPC
# ---------------------------------------------------
# $Rev: 212 $
# $LastChangedDate: 2011-01-24 12:08:32 +0100 (Mon, 24 Jan 2011) $
# $LastChangedBy: mpay $ 
#
# To compile you must load these modules:
#   module load pgi/2013
#   module load netcdf/pgi-4.2
#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

# compiler settings   
export FC        = pgf90
export FCFLAGS   =  -I$(IBMLIB_DIR) 
export FPPFLAGS  = 
export NETCDF_DIR = /appl/netcdf/pgi/4.2

export CC        = gcc    # why not pgcc ??
export AR        = ar
export RANLIB    = ranlib


MODDIRS = $(NETCDF_DIR)/include $(PHYSICAL_FIELDS_DIR) $(PARTICLE_STATE_DIR) $(TASK_DIR) $(OUTPUT_WRITER_DIRS)
FCFLAGS += $(addprefix -I,$(MODDIRS)) 

# linker settings
LINKFLAGS = -i4
LINKLIBS  += -L$(NETCDF_DIR)/lib -lnetcdff -lnetcdf -lcurl 


export LITTLE_ENDIAN = -convert little_endian
export BIG_ENDIAN    = -convert big_endian

