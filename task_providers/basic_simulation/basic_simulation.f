ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Standard particle tracking simulation 
c     ---------------------------------------------------
c
c     An almost minimal template
c     
c     setup and run an ensemble defined in input file 
c     task_providers/basic_simulation/basic_simulation_test_input.txt
c     writes trajectory of particle 1 to file trajectory1
c     writes final positions of ensemble to file finalpos 
c
c     to run it type: 
c        ibmrun task_providers/basic_simulation/basic_simulation_test_input.txt
c    
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
      use particle_state
      implicit none
c     ------------ declarations ------------      
      integer      :: idum4(4)
      type(clock),target  :: start_time, end_time
      type(clock),pointer :: current_time
      type(particle_ensemble)     :: par_ens
      type(emission_box), pointer :: emitboxes(:)
      integer :: year, month, day, julday,last,second_of_day
      integer :: i,ix,iy,iz,dum(79),ilast,iunit, isec
      integer :: idum(4), timeout, istep, ntracers
      real    :: time_step, now,xyz(3)
      real    :: amp,xcen,ycen,vdiff,hdiff
      character*256 :: filename
      type(clock) :: aclock

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
      open(44,file="trajectory1")

      open(72,file="Depth_particles")

      open(78,file="Particules_properties")
      write(78,*) "age, nb alage"

      do while (compare_clocks(current_time, end_time)
     +                         *nint(time_step) <= 0)             ! opposite for forward/backward simulation
         write(*,372) istep
         call update_physical_fields()
c        -------- propagate tracers  --------        
         call generate_particles(par_ens, emitboxes, time_step)
         call update_particles(par_ens, time_step)
c        do i = 1, 10, 1
c         write(78,*) istep, par_ens%state_stack(:)
c         write(78,*) istep, par_ens%state_stack
c        end do
c        -------- write tracer state --------
c        -------- loop control       --------
         call add_seconds_to_clock(current_time, nint(time_step))
         istep = istep + 1
         call get_last_particle_number(par_ens, last)
c         call set_clock_4i(aclock, year, month, day, second_of_day)
         call get_julian_day(current_time, julday)
         call get_second_in_day(current_time, isec)
         if (mod(julday,7)==0 .and. isec == 1) then
           do i=1,last
            write(78,*) istep, julday, i, i, i
            call write_state_attributes(par_ens%state_stack(i))
            call get_particle_position(get_particle(par_ens,i),xyz)
            write(72,*) istep,julday, i, xyz
           end do
         end if
         if (last>0) then
            call get_particle_position(get_particle(par_ens,1),xyz)
            write(44,*) istep, xyz
         endif
c        ---- quick way to dump e.g. vertical position of particles 1-10 for file fort.71        
c         if (last>9) then                ! only write, if particles 1-10 are actually emitted 
c            write(71,*) istep*time_step/60.0,
c     +                  par_ens%space_stack(1:10)%position(3)
c         endif
      enddo
      write(*,*) last
c     -------- dump final simulation output --------   
c      call write_ensemble(par_ens)
      open(45,file="finalpos")
      do i=1,last
         call get_particle_position(get_particle(par_ens,i),xyz)
         write(45,*) xyz
      enddo
      
      close(44)
      close(45)
 372  format(25("*"), " main time loop step ",i5.5, " ", 25("*"))
 377  format(f6.2, 1x, 3f8.3, 1x, i2, 1x, f8.3,i3)
c     ----------------- close down ---------------------------
     
      call close_physical_fields()
      call close_particles()

      write(*,*) "normal termination of simulation"

      end program
