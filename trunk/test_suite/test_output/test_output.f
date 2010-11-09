      program tracker
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Standard particle tracking simulation
c
c     make tracker
c     tracker task_providers/basic_simulation_test_input.txt
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use input_parser
      use time_tools
      use run_context
      use physical_fields
      use particles
      use netcdf_writer
      use ascii_writer
      use output
      use netcdf
    
      implicit none
c     ------------ declarations ------------      
      integer      :: idum4(4)
      type(clock),target  :: start_time, end_time
      type(clock),pointer :: current_time
      type(particle_ensemble)     :: par_ens
      type(emission_box), pointer :: emitboxes(:)
      integer :: year, month, day, julday,last
      integer :: i,ix,iy,iz,dum(79),ilast,iunit
      integer :: idum(4), timeout, istep, old_last
      real    :: time_step, now,xyz(3)
      real    :: amp,xcen,ycen,vdiff,hdiff
      character*256 :: filename

      type(ascii_output_file) :: birth_certs
      type(output_var),target,allocatable :: out_vars(:)
      type(netcdf_output_file) :: ncdf_trajs

      type(timer) :: start_timer
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
      call init_particles()
      call update_physical_fields()  ! need topology to verify emission boxes

      call setup_ensemble_from_file(par_ens, "emitbox", emitboxes)
      
c     ------------   setup output files  ------------   
      allocate(out_vars(7))   
      call construct(out_vars(1),"lat",var_type=NF90_SHORT,
     +      range=(/114.0,115.0/))
      call construct(out_vars(2),"lon",var_type=NF90_SHORT,  !Override default formatting
     +      range=(/0.0,2.0/))
      call construct(out_vars(3),"POSIX",var_type=NF90_INT)
      call construct(out_vars(4),"lat",fmt="(f7.4)")   !Override default formatting
      call construct(out_vars(5),"lon",fmt="(f7.4)")  !Override default formatting
      call construct(out_vars(6),"tracerID",var_type=NF90_INT)
      call construct(out_vars(7),"depth")
      birth_certs = ascii_output_file("bcs.txt",
     +    out_vars(3:7),char(9))   !Birth certificates
      ncdf_trajs = netcdf_output_file("trajs.nc",out_vars(3),
     +     get_ensemble_size(par_ens),
     +    out_vars(6),out_vars(1:2)) !Netcdf trajectory
      call init_output(birth_certs)
      call init_output(ncdf_trajs)

      current_time => get_master_clock()
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
            write(*,*) "Number of particles = ", last-old_last     
            call tic(start_timer)
!               call write_particles(birth_certs,par_ens,
!     +             (/(I,I=old_last+1,last)/))           !Implied do loop from old_last+1 to last
            call toc(start_timer,"Write particles")
         endif
c        -------- propagate tracers  --------
         call update_particles(par_ens, time_step)
c        -------- write tracer state -------- 
         call tic(start_timer)
         call write_frame(ncdf_trajs,par_ens)  
         call toc(start_timer,"Write frame")
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
