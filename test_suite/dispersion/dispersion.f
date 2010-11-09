      program tracker
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Standard particle tracking simulation
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use input_parser
      use time_tools
      use run_context
      use physical_fields
      use particles
     
      implicit none
c     ------------ declarations ------------      
      type(clock),target  :: release_time, end_time
      type(clock),pointer :: current_time
      type(particle_ensemble)     :: par_ens
      type(emission_box), pointer :: emitboxes(:)
      integer      :: idum4(4), istep,last,i, ipar,iunit,dir,nsteps
      real         :: xyz(3),time_step
      character*999 :: output_files(2)

c     ------------   show time starts  ------------
      call init_run_context()
      
c     ------------ initialise system --------------
      call random_seed()
      call init_physical_fields()
      call init_particles()
      call update_physical_fields()  ! need topology to verify emission boxes

c     ------------ read parameters from file --------------
      call read_control_data(simulation_file, "particle_time_step", 
     +                       time_step)
      call read_control_data(simulation_file,"backward_filename",
     +    output_files(1))
      call read_control_data(simulation_file,"forward_filename",
     +    output_files(2))
      call read_control_data(simulation_file, "emitbox", idum4)
      call set_clock(release_time,idum4(1),idum4(2),idum4(3),idum4(4))
      call read_control_data(simulation_file,"nsteps",nsteps)
     
c     ------------ loop forwards, then backwards --------------
      do dir=1,2

c     Setup ensemble      
      call setup_ensemble_from_file(par_ens, "emitbox", emitboxes)

c     Setup output file
      call find_free_io_unit(iunit)
      open(unit=iunit,file=trim(output_files(dir)))
      
c     Set end time and run direction
      time_step = abs(time_step)*(-1)**dir
      call set_master_clock(release_time)
      call set_clock(end_time,release_time) 
      call add_seconds_to_clock(end_time,nint(time_step*nsteps))
      current_time => get_master_clock()
      
c     =====================  main time loop =====================
      istep   = 0
      do while (compare_clocks(current_time, end_time)*time_step <= 0)
         write(*,372) dir, istep
         write(*,*) get_datetime(current_time)
         call update_physical_fields()
         
c        -------- release new tracers if needed  --------
         call generate_particles(par_ens, emitboxes, time_step)
         call get_last_particle_number(par_ens, last) 

c         -------- dump particle positions --------
         do i=1,last
            call get_particle_position(get_particle(par_ens,i),xyz)
            write(iunit,*) istep*time_step, xyz
         enddo  

c        -------- propagate tracers  --------
         call update_particles(par_ens, time_step)


c        -------- loop control       --------
         call add_seconds_to_clock(current_time, nint(time_step))
         istep = istep + 1
      enddo
      
      close(iunit)
      enddo
      
 372  format(25("*"), "loop ",i2, ", step ",i5.5, " ", 25("*"))
 
c     ----------------- close down ---------------------------
     
      call close_physical_fields()
      call close_particles()
      
      write(*,*) "normal end of simulation"

      end program
