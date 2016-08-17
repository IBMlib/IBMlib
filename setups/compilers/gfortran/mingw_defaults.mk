#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
# ----------------------------------------------------------------
# MinGW cross compiler (Linux -> Windows) default flags 
# 
# Currently distributed at https://sourceforge.net/projects/mingw
# Project home page: http://www.mingw.org/  (partially outdated)
# ----------------------------------------------------------------
# $Rev: 212 $
# $LastChangedDate: 2011-01-24 12:08:32 +0100 (Mon, 24 Jan 2011) $
# $LastChangedBy: mpay $ 
#
# This example illustrates how to include/link to a custom NetCDF build 
# export directive faciliates that variables are passed to sub-make
#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

# compiler settings   
export FC        =  /usr/bin/x86_64-w64-mingw32-gfortran     # 32 bit host to 64 bit target
export FCFLAGS   = -std=gnu -fall-intrinsics -I$(IBMLIB_DIR) # Fortran90 as standard can not be selected, -i4 is standard
export FPPFLAGS  = -fall-intrinsics                          # preprocessing automaticalyy invoked for .f, otherwise apply -cpp 
export CC        = /usr/bin/x86_64-w64-mingw32-gcc
export AR        = /usr/bin/x86_64-w64-mingw32-ar
export RANLIB    = /usr/bin/x86_64-w64-mingw32-ranlib

export NETCDF    = /home/asbjorn/DTU/Ballastvand_SRA/IBMlib_port_to_windows/mingw_netcdf/working/NETCDF

MODDIRS = $(PHYSICAL_FIELDS_DIR) $(PARTICLE_STATE_DIR) $(TASK_DIR) $(OUTPUT_WRITER_DIRS) $(NETCDF)/include
FCFLAGS += $(addprefix -I,$(MODDIRS))  

# linker settings
LINKFLAGS =    -static -fall-intrinsics  # -i4 is standard
LINKLIBS  +=    -L$(NETCDF)/lib -lnetcdff -lnetcdf 


export LITTLE_ENDIAN = -fconvert=little-endian
export BIG_ENDIAN    = -fconvert=big-endian
