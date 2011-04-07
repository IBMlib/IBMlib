##################################################################
#  ---------------------------------------------------
#  Makefile for running IBMlib tests
#  ---------------------------------------------------
#  $Rev$
#  $LastChangedDate$
#  $LastChangedBy$ 
#  ---------------------------------------------------
#  Notes: 
#    Makefile must appear as in .PHONY for reasons beyond comprehension
#    make with debugging info: make -d TARGET
#
#############################################################
.PHONY: fulltest FORCE clean variables $(FULLTEST_SET) Makefile

BASH         = bash
FULLTEST_SET = minimal_simulation test_output dispersion vertical_homogenity \
               interface_tests/particle interface_tests/pbi


fulltest: $(FULLTEST_SET)
	@echo "##################  result of FULLTEST ##################"
	@cat $(addsuffix /test_summary, $^)
	@echo "##################   end of FULLTEST   ##################"

clean: 
	@echo cleaning up
	@rm -f $(addsuffix /test_summary, $(FULLTEST_SET))
	@rm -f $(addsuffix /*.o, $(FULLTEST_SET))
	@rm -f $(addsuffix /*.a, $(FULLTEST_SET))
	@rm -f $(addsuffix /*.mod, $(FULLTEST_SET))


#
#  Apply default match-anything rule for performing a test ---
#  identified by the corresponding directory name. Test can be invoked independently
#  by typing "make <name_of_test>"
#
%: FORCE
	@echo "\n>>> performing test" $@ "<<<\n"
	@$(BASH) -c "\[  -d $@ \] && cd $@; runtest || echo 'Directory $@ not found'"


variables:
	@clear
        @echo "(FULLTEST_SET = " $(FULLTEST_SET)

FORCE: