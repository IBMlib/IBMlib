##################################################################
#  ---------------------------------------------------
#  Makefile for running IBMlib tests
#  ---------------------------------------------------
#  $Rev$
#  $LastChangedDate$
#  $LastChangedBy$ 
#  ---------------------------------------------------
#
#  Invoking test sets:
#
#    1) It is recommended to apply Setups/Compilers/FC/compiler_defaults_safe.mk
#       as $(IBMLIB_DIR)/compiler_defaults.mk
#    2) In $(IBMLIB_DIR) invoke: make fulltest
#    3) In this directory invoke: make fulltest or any TESTNAME
#   
#
#  Prerequisites:
#    It is assumed the standard Linux scripting languages are available wtihin 
#    the binary search path as well, in addition to prerequisites to IBMlib
#
#
#  Protocol for new tests:
#
#    1) A test is associated with a directory TESTNAME having the name of the test
#
#    2) The test must be executed by an (executable) script TESTNAME/runtest.
#       The script must be any standard Linux scripting language (*sh,python,perl)
#       The test invoked by TESTNAME/runtest must be self contained and able to
#       be invoked independtly. All data and sources must be available in TESTNAME 
#       or part of the resident IBMlib distribution. 
#
#    3) The test script must summarize the reults of the test in TESTNAME/test_summary
#       The test_summary should be short and conclusive and end with either line
#
#       "result TESTNAME: test OK"      or 
#       "result TESTNAME: test FAILED" (... or just crash)
#
#       In other words TESTNAME/runtest TESTNAME/test_summary are mandatory and reserved
#  
#    4) The code should preferably build into the executable ibmrun_test
#       to avoid accidental confusion with the normal name ibmrun for the
#       executable object
#
#  Other parts of test_suite:
#    * test_pbi is a stand alone test usefull for testing a new PBI
#
#
#  Notes: 
#    Makefile must appear as in .PHONY for reasons beyond comprehension
#    make with debugging info: make -d TARGET
#
#  TODO: 
#   minimal_simulation:     OK
#   test_output:            OK
#   dispersion               
#   interface_tests/particle 	       
#   interface_tests/pbi 
#   
#
#############################################################
.PHONY: fulltest FORCE clean variables $(FULLTEST_SET) Makefile

BASH         = bash
FULLTEST_SET = minimal_simulation test_output dispersion 
#              interface_tests/particle \
#	       interface_tests/pbi


fulltest: $(FULLTEST_SET)
	@echo "##################  result of FULLTEST ##################"
	@cat $(addsuffix /test_summary, $^)

clean: 
	@echo cleaning up
	@rm -f $(addsuffix /test_summary, $(FULLTEST_SET))

# --- Apply default match-anything rule for performing a test ---
%: FORCE
	@echo "\n>>> performing test" $@ "<<<\n"
	@$(BASH) -c "\[  -d $@ \] && cd $@; runtest || echo 'Directory $@ not found'"


variables:
	@clear
        @echo "(FULLTEST_SET = " $(FULLTEST_SET)

FORCE: