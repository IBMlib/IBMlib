      program tracker
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Test_output
c
c     Demonstrates how the output modules can be used to just write the 
c     final state of the ensemble. Particle properties written to file
c     is controlled from the input file by tags "final_state_ascii_data" (additively, in order).
c     Output file name is controleld by tag "final_state_ascii_file".
c
c     This task provider will inquire particle properties that are 
c     provided by e.g. biology_providers/passive
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use input_parser
      use time_tools
      use run_context
      use physical_fields
      use particles
      use ascii_writer
      use output
    
      implicit none
c     ------------ declarations ------------      
      type(clock),target  :: start_time, end_time
      type(clock),pointer :: current_time
      type(particle_ensemble)     :: par_ens
      type(emission_box), pointer :: emitboxes(:)
      type(ascii_output_file) :: final
     
      character*999 , allocatable :: out_vars(:)
      type(variable), allocatable :: vars(:)
      integer :: i, idum4(4), istep, last,old_last
      real    :: time_step

c     ------------   show time starts  ------------
      call init_run_context()

c     ------------   set clocks  ------------   
      call read_control_data(simulation_file, "start_time", idum4)
      call set_clock(start_time, idum4(1), idum4(2), idum4(3), idum4(4))
      call read_control_data(simulation_file, "end_time", idum4)
      call set_clock(end_time, idum4(1), idum4(2), idum4(3), idum4(4))
      call read_control_data(simulation_file, "particle_time_step", 
     +                       time_step)

c     ------------   init system  ------------   
      call init_physical_fields(start_time)
      current_time => get_master_clock()
      call init_particles()
      call update_physical_fields()  ! need topology to verify emission boxes
      call setup_ensemble_from_file(par_ens, "emitbox", emitboxes)
      
c     ------------   setup output files  ------------   

      call setup_output_from_file(final,simulation_file,
     +     "final_state_ascii_file", "final_state_ascii_data")
      

c     =====================  main time loop =====================
      istep   = 0
      do while (compare_clocks(current_time, end_time) <= 0)
         write(*,372) istep
         call update_physical_fields()
c        -------- release new tracers and record new tracers  --------
         call generate_particles(par_ens, emitboxes, time_step)         
c        -------- propagate tracers  --------
         call update_particles(par_ens, time_step)
c        -------- write trajectories --------
c        -------- loop control       --------
         call add_seconds_to_clock(current_time, nint(time_step))
         istep = istep + 1
      enddo
c     -------- write final frame --------
      call get_last_particle_number(par_ens, last) 
      write(*,*) "writing final state for ",last,"particles"
      call write_frame(final,par_ens)   
      
 372  format(25("*"), " main time loop step ",i5.5, " ", 25("*"))

c     ----------------- close down ---------------------------
      call close_physical_fields()
      call close_particles()

      write(*,*) "normal end of simulation"

      end program
