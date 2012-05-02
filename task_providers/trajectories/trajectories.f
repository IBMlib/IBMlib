ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Trajectories Generator
c     ---------------------------------------------------
c     $Rev$
c     $LastChangedDate$
c     $LastChangedBy$ 
c
c     Generic task provider to generate a set of trajectories. All 
c     configuration, including the selection of the output variables
c     and their frequency should take place through the configuration
c     file. System is intended to be generic and work for both forward
c     and backwards tracking. Use of multiple particle ensembles are
c     supported, allowing for the optimisation of storage
c     Task provider allows for the following types of outputs
c      <release_points>  ASCII file listing the release state of each 
c                        particle
c      <trajectories>    NetCDF file listing the value of key 
c                        variables as a function of time
c      <final_state>     ASCII file listing the value of chosen
c                        variables at completion of simulation
c
c     TODO: 
c     1. Multiple particle ensembles should write to separate files
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program tracker
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
      type(timer)         :: stopwatch
      type(particle_ensemble)     :: par_ens
      type(emission_box), pointer :: emitboxes(:)
      type(ascii_output_file)  :: final_state,birth_certs
      type(netcdf_output_file) :: trajs
      integer :: i, idum4(4), istep, last,old_last,seed_size
      integer, allocatable :: seed(:)
      integer :: opt_freq,i8(8)
      real    :: time_step
      logical :: write_traj,write_release,write_final
      logical :: write_traj_rng
      type(time_period) :: write_period
      character*200 :: ver_str 

      call tic(stopwatch)
c     ------------   Display version numbers  ------------
      write(*,*) "Trajectories task provider:  $Rev$"
      write(*,*) trim(get_pbi_version())
      write(*,*) trim(get_particle_version())

c     ------------   Initialise config file   ------------
      call init_run_context()

c     ------------   set clocks  ------------   
      call read_control_data(simulation_file, "start_time", idum4)
      call set_clock(start_time, idum4(1), idum4(2), idum4(3), idum4(4))
      call read_control_data(simulation_file, "end_time", idum4)
      call set_clock(end_time, idum4(1), idum4(2), idum4(3), idum4(4))
      call read_control_data(simulation_file, "time_step", 
     +                       time_step)

c     ------------   init system  ------------   
      call init_physical_fields(start_time)
      current_time => get_master_clock()
      call init_particles()

c     ------------ seed random number generator  ------------
      call random_seed(size=seed_size)    !get seed size
      allocate(seed(seed_size))
      if(count_tags(simulation_file,"random_seed")>0) then 
           call read_control_data(simulation_file, "random_seed", seed)
           call random_seed(put=seed)
      else
           call random_seed()
      endif
      call random_seed(get=seed)
      write (*,*) "random_seed = ", seed

c     ------------   setup particle ensemble(s) ------------   
      call setup_ensemble_from_file(par_ens, "emitbox", emitboxes)
      call update_physical_fields()  

c     ------------   setup output files  ------------   
      write_traj = count_tags(simulation_file,"trajectory_fname")/=0
      write_traj_rng=count_tags(simulation_file,"traj_output_range")/=0
      write_release = count_tags(simulation_file,"release_dat_fname")/=0
      write_final = count_tags(simulation_file,"final_state_fname")/=0
      if(write_release) then
         call setup_output_from_file(birth_certs,simulation_file,
     +      "release_dat_fname","release_var")
      endif
      if(write_final) then
         call setup_output_from_file(final_state,simulation_file,
     +      "final_state_fname","final_state_var")
      endif
      if(write_traj) then
         call read_control_data(simulation_file, 
     +                  "traj_output_frequency",opt_freq)
         call setup_output_from_file(trajs,par_ens,simulation_file,
     +      "trajectory_fname","traj_var")
      endif
      if(write_traj_rng) then
         call read_control_data(simulation_file, 
     +                  "traj_output_range",i8)
         call set_period(write_period,i8)
         write(*,*) "Writing trajectory data during period:"
         call write_time_period(write_period)
      endif

c     =====================  main time loop =====================
      istep   = 0
      do while (compare_clocks(current_time, end_time)
     *             *nint(time_step) <= 0)
         write(*,372)istep, get_datetime(current_time)
         call update_physical_fields()
c        -------- release new tracers and record new tracers  --------
         call get_last_particle_number(par_ens, old_last) 
         call generate_particles(par_ens, emitboxes, time_step)
         if(write_release) then
            call get_last_particle_number(par_ens, last) 
            if(last/=old_last) then !new tracers released!
            call write_particles(birth_certs,par_ens,
     +             (/(I,I=old_last+1,last)/))           !Implied do loop from old_last+1 to last
            endif
         endif
c        -------- propagate tracers  --------
         call update_particles(par_ens, time_step)
c        -------- write trajectories -------- 
         if(write_traj)  then
            if(MOD(istep,opt_freq)==0)  then
            if(.not. write_traj_rng) then
               call write_frame(trajs,par_ens)  
            else if(time_inside_period(current_time,write_period)) then
               call write_frame(trajs,par_ens)  
            endif
            endif
         endif
         
c        -------- loop control       --------
         call add_seconds_to_clock(current_time, nint(time_step))
         istep = istep + 1
      enddo
 372  format(25("*"), " step ",i6, ": ",a20," ", 25("*"))

c     ----------------- Write final state --------------------
      if(write_final) then
         call write_frame(final_state,par_ens)
      endif

c     ----------------- close down ---------------------------
      call close_physical_fields()
      call close_particles()

      call toc(stopwatch,"Simulation run time")
      write(*,*) "=========normal end of simulation========="


      end program
