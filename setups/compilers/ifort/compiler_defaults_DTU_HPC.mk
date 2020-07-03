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
#
# HPC configuration: 
#   module load intel/2017.4.196
#   module load netcdf-fortran/4.4.4-intel-2017update4
#
#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

# compiler settings   
export FC        = ifort
export FCFLAGS   = -i4  -I$(IBMLIB_DIR) # -e90 -i4 -error_limit 3 -I$(IBMLIB_DIR)
export FPPFLAGS  = -fpp 
export CC        = icc
export AR        = ar
export RANLIB    = ranlib

MODDIRS = $(PHYSICAL_FIELDS_DIR) $(PARTICLE_STATE_DIR) $(TASK_DIR) $(OUTPUT_WRITER_DIRS)
FCFLAGS += $(addprefix -I,$(MODDIRS))  
FCFLAGS += -I/usr/local/include  -I/appl/netcdf-fortran/4.4.4/intel-2017update4/include # should contain netcdf.mod


# linker settings
LINKFLAGS = -i4
LINKLIBS  += -L/usr/local/lib  -L/appl/netcdf-fortran/4.4.4/intel-2017update4/lib   # 
LINKLIBS  += -lnetcdff -lnetcdf  -lhdf5_hl -lhdf5 -lz  # netcdf

export LITTLE_ENDIAN = -convert little_endian
export BIG_ENDIAN    = -convert big_endian
