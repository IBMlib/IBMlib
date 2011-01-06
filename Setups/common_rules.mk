#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#     ---------------------------------------------------
#     Makefiles common rules
#     ---------------------------------------------------
#     $Rev: 147 $
#     $LastChangedDate: 2010-11-19 00:33:57 +0100 (Fri, 19 Nov 2010) $
#     $LastChangedBy: mpay $ 
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
