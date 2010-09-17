##################################################################
#
#               Basic makefile for IBMlib (requires gmake) 
#         
#  This makefile builds IBMLIB_CORE and delegates out the update  
#  of modules PHYSICAL_FIELDS/PARTICLE_STATE/TASK and makes to 
#  final linking to all objects into the final executable EXECUTABLE.
#  IBMLIB has six levels giving the appropriate build order (and allowed
#  directions of use association)
#
#   1) TASK 
#   2) particles.mod
#   3) PARTICLE_STATE 
#   4) particle_tracking.mod
#   5) PHYSICAL_FIELDS 
#   6) IBMLIB_BASE (incl external included tools)     
#
#  Each actual module PHYSICAL_FIELDS/PARTICLE_STATE/TASK must  be associated
#  with a separate directory containing one makefile providing the (minimal) targets
#  listed below.
#  One level is allowed to use associate to the same or a lower level
#  One level is not allowed to use associate to a higher level
#  One level is only allowed to build objects at its own level (distributed makefiles)
#  This means that .mod files from other levels should be considered external
#  like system installations. This means that makefile for PHYSICAL_FIELDS/PARTICLE_STATE/TASK
#  should suppress updates of external components and that root makefile
#  is responsible for overall build order 6->5->4->3->2->1 above. 
#  Each makefile of  may PHYSICAL_FIELDS/PARTICLE_STATE/TASK - or may not - include common.mk
#
#  Update responsibilities:
#  	$(IBMLIB_DIR)/Makefile:          IBMLIB_OBJECTS   EXECUTABLE 
#       $(PHYSICAL_FIELDS_DIR)/Makefile: PHYSICAL_FIELDS
#       $(PARTICLE_STATE_DIR)/Makefile:  PARTICLE_STATE
#       $(TASK_DIR)/Makefile:            TASK
#  ---------------------------------------------------------------------------
#  module PHYSICAL_FIELDS
#     rooted in PHYSICAL_FIELDS_DIR. In this directory, there
#     should be a makefile updating the targets:
#        physical_fields.mod (F90 module interface, in directory PHYSICAL_FIELDS_DIR)
#        physical_fields.a   (all compiled objects of module, in directory PHYSICAL_FIELDS_DIR)
#        clean
#     An (optional) makefile link_opt.mk in PHYSICAL_FIELDS_DIR may define the following 
#     variables for link options to be used for the final stage linking:
#        LINKFLAGS_PHYSICAL
#        LINKLIBS_PHYSICAL 
#  ---------------------------------------------------------------------------   
#  module PARTICLE_STATE
#     rooted in PARTICLE_STATE_DIR. In this directory, there
#     should be a makefile updating the targets:
#        particle_state.mod (F90 module interface, in directory PARTICLE_STATE_DIR)
#        particle_state.a   (all compiled objects of module, in directory PARTICLE_STATE_DIR)
#        clean    
#     An (optional) makefile link_opt.mk in PHYSICAL_FIELDS_DIR may define the following 
#     variables for link options to be used for the final stage linking:
#        LINKFLAGS_STATE
#        LINKLIBS_STATE      
#  ---------------------------------------------------------------------------   
#  module TASK. 
#     rooted in TASK_DIR  In this directory, there
#     should be a makefile updating the targets:
#        task.a (all compiled objects INCLUDING the main program, in directory PARTICLE_STATE_DIR)
#        clean
#     An (optional) makefile link_opt.mk in PHYSICAL_FIELDS_DIR may define the following 
#     variables for link options to be used for the final stage linking:
#        LINKFLAGS_TASK
#        LINKLIBS_TASK     
#  ---------------------------------------------------------------------------
#  Restrictions:
#     1) To have an intelligeble Makefile these rule apply for code organisation
#        so that a .mod.f rule can be applied for updating F90 modules:
# 
#        1a) one-module-one-file: only one module in each file
#        1b) head names of file and module match, i.e. yyy.f expresses module yyy
#        1c) a simple consistent case convention for module names 
#            (so module_dependence_scan.py can express dependency rules)
#        
#  TODO: enforce strict directions of use association
#        online/offline configuration of task
#
#  About order of objects at linking: "The traditional behavior of linkers is to search 
#  for external functions from left to right in the libraries specified on the command line. 
#  This means that a library containing the definition of a function should appear after any 
#  source files or object files which use it ...When several libraries are being used, 
#  the same convention should be followed for the libraries themselves."
#
#  ifort note: apparently this problem with ifort can just be handled
#  by duplicating link objects once ...


##################################################################
#              M a i n   c o n f i g u r a t i o n               #
##################################################################
EXECUTABLE          = ibmrun   # name of the executable to build
PHYSICAL_FIELDS_DIR = $(IBMLIB_DIR)/oceanography_providers/vortex1
PARTICLE_STATE_DIR  = $(IBMLIB_DIR)/biology_providers/passive
TASK_DIR            = $(IBMLIB_DIR)/task_providers/test
##################################################################

export IBMLIB_DIR  = $(shell pwd)


IBMLIB_EXTERNAL =  libtime/libtime77.a
IBMLIB_BASEOBJ  =  constants.o  grid_interpolations.o  input_parser.o  random_numbers.o  runtime_tools.o\
                   time_tools.o geometry.o   run_context.o  string_tools.o  $(IBMLIB_EXTERNAL)
IBMLIB_BASEMOD  =  constants.mod  input_parser.mod  random_numbers.o  time_tools.mod    run_context.mod 
IBMLIB_OBJECTS  = $(IBMLIB_BASEOBJ)  physical_fields.a  particle_tracking.o particle_state.a particles.o  task.a 


#     Main task of this Makefile: EXECUTABLE - should be first target, to appear as default
# --- currently TASK appears as task.a - at some point we may homogenize 
#     the setup. Dependencies order reflect strict build order of IBMlib
#
#     About order of objects at linking: the  linker  searches  and processes libraries and 
#     object files in the order they are specified, you should specify 
#     providing functions AFTER the last object file it applies. This is 
#     exactly opposite the make order. 
#     Therefore hack: dublicate IBMLIB_OBJECTS to ensure a provider is also after using function ...

$(EXECUTABLE): $(IBMLIB_OBJECTS)
	$(FC)  $(IBMLIB_OBJECTS) $(IBMLIB_OBJECTS) $(LINKFLAGS) $(LINKLIBS) -o $(EXECUTABLE)


include $(IBMLIB_DIR)/common.mk   # shared information between makefiles

# 
# --- collect additional optional linkflags / link options from PHYSICAL_FIELDS/PARTICLE_STATE/TASK 
#

LINKFLAGS = -i4 
LINKLIBS  = 
-include $(PHYSICAL_FIELDS_DIR)/link_opt.mk   # optional include 
-include $(PARTICLE_STATE_DIR)/link_opt.mk    # optional include 
-include $(TASK_DIR)/link_opt.mk              # optional include 
LINKFLAGS += $(LINKFLAGS_PHYSICAL)  $(LINKFLAGS_STATE)  $(LINKFLAGS_TASK)
LINKLIBS  += $(LINKLIBS_PHYSICAL)   $(LINKLIBS_STATE)   $(LINKLIBS_TASK)

#
# --- and here below goes other targets of this makefile ---
#

physical_fields.a: FORCE  $(IBMLIB_BASEMOD) 
	cd $(PHYSICAL_FIELDS_DIR); make physical_fields.mod physical_fields.a
	ln -sf $(PHYSICAL_FIELDS_DIR)/physical_fields.a
	ln -sf $(PHYSICAL_FIELDS_DIR)/physical_fields.mod

particle_state.a:  FORCE  $(IBMLIB_BASEMOD)  physical_fields.a  particle_tracking.o
	cd $(PARTICLE_STATE_DIR);  make particle_state.mod  particle_state.a
	ln -sf $(PARTICLE_STATE_DIR)/particle_state.a
	ln -sf $(PARTICLE_STATE_DIR)/particle_state.mod

task.a: FORCE  $(IBMLIB_BASEMOD)  physical_fields.a  particle_tracking.mod particle_state.a particles.mod 
	cd $(TASK_DIR);  make task.a 
	ln -sf $(TASK_DIR)/task.a 


# take libtime/Makefile diagnostic on whether libtime77 has been untarred, hmm...


libtime/libtime77.a: FORCE  libtime/Makefile  
	cd libtime;  make libtime77.a

libtime/Makefile:
	cd libtime; tar xvfz libtime.tar.gz; rm -f Makefile; \
        ln -s Makefile_adapted_asc Makefile


clean: FORCE
	/bin/rm -f *.o *.a *.mod dependences.mk $(EXECUTABLE)
	cd libtime; make cleanall
	cd $(PHYSICAL_FIELDS_DIR); make clean
	cd $(PARTICLE_STATE_DIR);  make clean
	cd $(TASK_DIR);            make clean


FORCE:


include $(IBMLIB_DIR)/common.mk   # gmake settings
