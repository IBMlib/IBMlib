ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Minimal particle tracking simulation test
c     ---------------------------------------------------
c
c     An almost minimal example - do nothing but    
c     setup and run an ensemble defined in input file 
c     task_providers/basic_simulation/basic_simulation_test_input.txt
c     
c     Test success criterium: the code builds and the simulation completes normally
c     ---------------------------------------------------
c     $Rev$
c     $LastChangedDate$
c     $LastChangedBy$ 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program tracker
      use input_parser
      use time_tools
      use run_context
      use physical_fields
      use particles
     
      implicit none
c     ------------ declarations ------------      
      integer      :: idum4(4)
      type(clock),target  :: start_time, end_time
      type(clock),pointer :: current_time
      type(particle_ensemble)     :: par_ens
      type(emission_box), pointer :: emitboxes(:)
      integer :: year, month, day, julday,last
      integer :: i,ix,iy,iz,dum(79),ilast,iunit
      integer :: idum(4), timeout, istep, ntracers
      real    :: time_step, now,xyz(3)
      real    :: amp,xcen,ycen,vdiff,hdiff
      character*256 :: filename
      
c     ------------   show time starts  ------------
      call init_run_context()

c.....set clocks 
      call read_control_data(simulation_file, "start_time", idum4)
      call set_clock(start_time, idum4(1), idum4(2), idum4(3), idum4(4))
      call read_control_data(simulation_file, "end_time", idum4)
      call set_clock(end_time, idum4(1), idum4(2), idum4(3), idum4(4))
      call read_control_data(simulation_file, "particle_time_step", 
     +                       time_step)

      call init_physical_fields(start_time)
      call init_particles()
      call update_physical_fields()  ! need topology to verify emission boxes

      call setup_ensemble_from_file(par_ens, "emitbox", emitboxes)
      
      current_time => get_master_clock()
c     =====================  main time loop =====================
      istep   = 0
      do while (compare_clocks(current_time, end_time) <= 0)
         call update_physical_fields()
c        -------- propagate tracers  --------        
         call generate_particles(par_ens, emitboxes, time_step)
         call update_particles(par_ens, time_step)
c        -------- write tracer state -------- 
c        -------- loop control       --------
         call add_seconds_to_clock(current_time, nint(time_step))
         istep = istep + 1
         call get_last_particle_number(par_ens, last) 
      enddo
c     -------- dump final simulation output --------   

c     ----------------- close down ---------------------------
      call close_physical_fields()
      call close_particles()

c     ------- normal simulation ends: write test report -------
      open(31,file="test_summary") 
      write(31,*) "result minimal_simulation: test OK"
      close(31)
c
      end program
