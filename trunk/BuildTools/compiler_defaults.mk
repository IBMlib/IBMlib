# ==========================================
# compiler settings   
# export directive faciliates that variables are passed to sub-make
# (unless locally overwritten)
# ==========================================

# compiler settings   
export FC        = ifort
export FCFLAGS   =  -e90 -i4 -error_limit 3 -I$(IBMLIB_DIR)  
export FPPFLAGS  = -fpp 

# linker settings
LINKFLAGS = -i4
LINKLIBS  += -L/usr/local/include 


