#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#     ---------------------------------------------------
#     Makefiles common rules
#     ---------------------------------------------------
#     $Rev$
#     $LastChangedDate$
#     $LastChangedBy$ 
#
#     Common rules shared by Makefiles
#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

# ============== implicit rules ==============
%.o : %.f
	$(FC) -c $(FCFLAGS) $(FPPFLAGS) $< -o $@

%.mod : %.f
	$(FC) -c $(FCFLAGS) $(FPPFLAGS) $< -o $*.o

# ============== module dependeces ==============
dependences.mk: *.f
	@$(BUILD_TOOLS)/module_dependence_scan.py $(EXTRES) *.f > dependences.mk

-include dependences.mk
