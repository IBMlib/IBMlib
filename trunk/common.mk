# ==========================================
# IBMLIB relative references
# ==========================================

export BUILD_TOOLS = $(IBMLIB_DIR)/BuildTools
export VPATH       = $(IBMLIB_DIR)  # make search path for src/obj    


# ==========================================
# compiler settings   
# export directive faciliates that variables are passed to sub-make
# (unless locally overwritten)
# ==========================================

export FC        = ifort
# big_endian needed to read DMI data
export FCFLAGS   = -g -e90 -i4 -error_limit 3 -convert big_endian   # -c -O2 -w 
FCFLAGS         += -I $(IBMLIB_DIR)
export FPPFLAGS  = -fpp 
export VPATH     = $(IBMLIB_DIR)  # make search path for src/obj    


# ==========================================
#           common gmake settings
# ==========================================
# 
.SUFFIXES: .mod   # add this suffix to set

# ============== implicit rules ==============
%.o : %.f
	$(FC) -c $(FCFLAGS) $(FPPFLAGS) $< -o $@

%.mod : %.f
	$(FC) -c $(FCFLAGS) $(FPPFLAGS) $< -o $*.o

# ============== module dependeces ==============
dependences.mk: *.f
	$(BUILD_TOOLS)/module_dependence_scan.py $(EXTRES) *.f > dependences.mk

include dependences.mk
