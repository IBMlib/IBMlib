#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
# ----------------------------------------------------------------
# gfortran (from GNU Compiler Collection)  default flags 
# ----------------------------------------------------------------
# $Rev: 212 $
# $LastChangedDate: 2011-01-24 12:08:32 +0100 (Mon, 24 Jan 2011) $
# $LastChangedBy: mpay $ 
#
# compiler settings   
# export directive faciliates that variables are passed to sub-make
# (unless locally overwritten)
#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

# compiler settings   
export FC        = gfortran
export FCFLAGS   = -std=gnu -fall-intrinsics -I$(IBMLIB_DIR) #  Fortran90 as standard can not be selected, -i4 is standard
export FPPFLAGS  = -fall-intrinsics                          # preprocessing automaticalyy invoked for .f, otherwise apply -cpp 

MODDIRS = $(PHYSICAL_FIELDS_DIR) $(PARTICLE_STATE_DIR) $(TASK_DIR) $(OUTPUT_WRITER_DIRS)
FCFLAGS += $(addprefix -I,$(MODDIRS))  

# linker settings
LINKFLAGS  = -fall-intrinsics   # -i4 is standard
LINKLIBS  +=                    # -L/usr/local/include 


export LITTLE_ENDIAN = -fconvert=little-endian
export BIG_ENDIAN    = -fconvert=big-endian
