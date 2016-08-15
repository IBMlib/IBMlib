#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
# ---------------------------------------------------
# Intel Fortran (ifort) Compiler default flags 
# ---------------------------------------------------
# $Rev: 212 $
# $LastChangedDate: 2011-01-24 12:08:32 +0100 (Mon, 24 Jan 2011) $
# $LastChangedBy: mpay $ 
#
# compiler settings   
# export directive faciliates that variables are passed to sub-make
# (unless locally overwritten)
#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

# compiler settings   
export FC        = ifort
export FCFLAGS   =  -e90 -i4 -error_limit 3 -I$(IBMLIB_DIR)
export FPPFLAGS  = -fpp 
export CC        = gcc
export AR        = ar


MODDIRS = $(PHYSICAL_FIELDS_DIR) $(PARTICLE_STATE_DIR) $(TASK_DIR) $(OUTPUT_WRITER_DIRS)
FCFLAGS += $(addprefix -I,$(MODDIRS))  
FCFLAGS += -I/usr/local/include   # should contain netcdf.mod

# linker settings
LINKFLAGS = -i4
LINKLIBS  += -L/usr/local/include 
LINKLIBS  += -lnetcdff -lnetcdf  -lhdf5_hl -lhdf5 -lz  # netcdf

export LITTLE_ENDIAN = -convert little_endian
export BIG_ENDIAN    = -convert big_endian
