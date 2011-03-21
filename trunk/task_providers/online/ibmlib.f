ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Module cast of standard minimal particle tracking simulation
c     for online demonstration usage
c
c     Unresolved issues: 
c
c        1) init_run_context assumes first runtime argument is an
c           input file for IBMlib ...
c
c        2) what fill value is used? (-> fill_value_real in online_pbi.f)
c
c        3) is it possible to transfer grid info generically ?           
c
c
c     $Rev:  $
c     $LastChangedDate: $
c     $LastChangedBy: $ 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module ibmlib
      use input_parser
      use time_tools
      use run_context
      use physical_fields
      use particles
      implicit none

c     --------  define scope --------
      private
      public :: physical_data_bucket                 ! reexport from physical_fields
      public :: transfer_data_to_ibmlib              ! reexport from physical_fields
      public :: init_ibmlib, ibm_step, close_ibmlib  ! defined in this module

c     --------  module data section  --------
     
      type(clock),target          :: start_time
      type(clock),pointer         :: current_time
      type(particle_ensemble)     :: par_ens
      type(emission_box), pointer :: emitboxes(:)
      integer                     :: istep    ! step counter for IBM_step
      integer                     :: trajfile
      integer                     :: finalfile
      logical                     :: ensemble_is_configured 

c     ----------------------------------
                  contains
c     ----------------------------------


      subroutine init_ibmlib()
c     -----------------------------------------------------
c     Initialization of simulation and necessary 
c     data structures
c
c     Defer initialization of the ensemble, since this 
c     generally requires that transfer_data_to_ibmlib and update_physical_fields
c     has been invoked (so that the current topography/hydrography
c     is available)
c     -----------------------------------------------------
      integer      :: idum4(4)
c     -------------------------------------------      
      call init_run_context()

c.....set clocks 
      call read_control_data(simulation_file, "start_time", idum4)
      call set_clock(start_time, idum4(1), idum4(2), idum4(3), idum4(4))
      call init_physical_fields(start_time)
      call init_particles()
      istep                  = 0
      ensemble_is_configured = .false.
      call find_free_IO_unit(trajfile)
      open(trajfile, file="trajectory1")
      call find_free_IO_unit(finalfile)
      open(finalfile, file="finalpos")
      end subroutine init_ibmlib
      



      subroutine ibm_step(time_step)
c     -----------------------------------------------------
c     Propagate IBM ensemble forward in time by time_step
c     Do not monitor simulation period, just propagate master_clock
c     by time_step
c     -----------------------------------------------------
      real, intent(in)  :: time_step

      real              :: xyz(3)
      integer           :: last
c     ------------------------------------------- 
      current_time => get_master_clock()      
      write(*,372) istep
 372  format(25("*"), " IBMlib main step ",i5.5, " ", 25("*"))      
      call update_physical_fields()
c     -------- check ensemble is ready  -------- 
      if (.not.ensemble_is_configured) then   
         call setup_ensemble_from_file(par_ens, "emitbox", emitboxes)
         ensemble_is_configured = .true.
      endif
c     -------- propagate tracers  --------        
      call generate_particles(par_ens, emitboxes, time_step)
      call update_particles(par_ens, time_step)
      call add_seconds_to_clock(current_time, nint(time_step))
      istep = istep + 1
c     -------- write tracer state time frame -------- 
      call get_last_particle_number(par_ens, last) 
      if (last>0) then
         call get_particle_position(get_particle(par_ens,1),xyz)
         write(trajfile,*) xyz
      endif
      
      end subroutine ibm_step





      subroutine close_ibmlib()
c     -----------------------------------------------------
c     Write final state and close down IBMlib
c     -----------------------------------------------------
      real              :: xyz(3)
      integer           :: i,last
c     -------------------------------------------
c      
c     ------  write final state  ------
c
      call get_last_particle_number(par_ens, last) 
      do i=1,last
         call get_particle_position(get_particle(par_ens,i),xyz)
         write(finalfile,*) xyz
      enddo
c     ------  close down         ------
      close(trajfile)
      close(finalfile)
      call close_physical_fields()
      call close_particles()
      write(*,*) "IBMlib: normal termination of simulation"
c
      end subroutine  close_ibmlib



      end module
