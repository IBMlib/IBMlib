      program tracker
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Test_output
c
c     Demonstrates how the output modules can be used to write the
c     results of a simulation. Two output files are generated
c     1. "release_data.txt" - An ASCII file containing the points in
c        time and space where each particle was released. This file is
c        generated using the "ascii_writer" module
c     2. "trajectories.nc" - A NetCDF file containing the lat and lot
c        of each particle at each time step. This file is generated
c        using the "netcdf_writer" module
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use input_parser
      use time_tools
      use run_context
      use physical_fields
      use particles
      use netcdf_writer
      use ascii_writer
      use output
    
      implicit none
c     ------------ declarations ------------      
      type(clock),target  :: start_time, end_time
      type(clock),pointer :: current_time
      type(particle_ensemble)     :: par_ens
      type(emission_box), pointer :: emitboxes(:)
      type(ascii_output_file) :: birth_certs
      type(netcdf_output_file) :: trajs
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
      allocate(out_vars(5))
      allocate(vars(size(out_vars)))
      out_vars(1)  = "tracerID"
      out_vars(2)  = "sourceBox"
      out_vars(3:4)=(/"lat","lon"/)
      out_vars(5)  = "depth"
      call get_metadata(out_vars,vars)
      call init_output(birth_certs,"release_data.txt",",",vars)
      call set_type(vars(3:4),"NF90_SHORT")
      call set_type(vars(5),"NF90_FLOAT")
      call set_var_range(vars(3),(/45.0,65.0/))
      call set_var_range(vars(4),(/-15.0,15.0/))
      call init_output(trajs,"trajectories.nc",
     +       vars(3:5),get_ensemble_size(par_ens))

c     =====================  main time loop =====================
      istep   = 0
      do while (compare_clocks(current_time, end_time) <= 0)
         write(*,372) istep
         call update_physical_fields()
c        -------- release new tracers and record new tracers  --------
         call get_last_particle_number(par_ens, old_last) 
         call generate_particles(par_ens, emitboxes, time_step)
         call get_last_particle_number(par_ens, last) 
         if(last/=old_last) then !new tracers released!
            call write_particles(birth_certs,par_ens,
     +             (/(I,I=old_last+1,last)/))           !Implied do loop from old_last+1 to last
         endif
c        -------- propagate tracers  --------
         call update_particles(par_ens, time_step)
c        -------- write trajectories -------- 
         call write_frame(trajs,par_ens)  
c        -------- loop control       --------
         call add_seconds_to_clock(current_time, nint(time_step))
         istep = istep + 1
      enddo
 372  format(25("*"), " main time loop step ",i5.5, " ", 25("*"))

c     ----------------- close down ---------------------------
      call close_physical_fields()
      call close_particles()

      write(*,*) "normal end of simulation"

      end program
