#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
# ---------------------------------------------------
# Intel Fortran (ifort) Compiler default flags 
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
export FCFLAGS   =  -e90 -i4 -error_limit 3 -I$(IBMLIB_DIR)  
export FPPFLAGS  = -fpp 

# linker settings
LINKFLAGS = -i4
LINKLIBS  += -L/usr/local/include 


