#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
# ---------------------------------------------------
# Intel Fortran (ifort) Compiler safe flags 
# ---------------------------------------------------
# $Rev: 147 $
# $LastChangedDate: 2010-11-19 00:33:57 +0100 (Fri, 19 Nov 2010) $
# $LastChangedBy: mpay $ 
#
# compiler settings   
# export directive faciliates that variables are passed to sub-make
# (unless locally overwritten)
#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

# compiler settings   
export FC        = ifort
export FCFLAGS   =  -e90 -i4 -error_limit 3 -I$(IBMLIB_DIR) -check bounds -check uninit -traceback  #Add bounds checking and tracebacks
export FPPFLAGS  = -fpp 

MODDIRS = $(PHYSICAL_FIELDS_DIR) $(PARTICLE_STATE_DIR) $(TASK_DIR) $(OUTPUT_WRITER_DIRS)
FCFLAGS += $(addprefix -I,$(MODDIRS))  


# linker settings
LINKFLAGS = -i4
LINKLIBS  += -L/usr/local/include 


export LITTLE_ENDIAN = -convert little_endian
export BIG_ENDIAN    = -convert big_endian
