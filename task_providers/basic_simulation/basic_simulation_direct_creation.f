ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Standard particle tracking simulation 
c     ---------------------------------------------------
c
c     setup and run an ensemble where particles are created by direct operator make_particles_direct
c     direct operator make_particles_direct
c     Input file example (only task part)
c     
c         task_providers/basic_simulation/basic_simulation_test_direct.txt
c
c     This example writes trajectory of particle 1 to file trajectory1
c     and writes final positions of ensemble to file finalpos 
c    
c     Using emission boxes, particles in a standard simulation are created by subroutine
c     setup_ensemble_from_file(...). Direct creation requireres invoking the ensemble
c     generator setup_ensemble(...) directly at start and create each particle by generate_particles
c     (generate_particles is overloaded with different argument sequences for alternative calls)
c    
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
      integer :: year, month, day, julday,last
      integer :: i,ix,iy,iz,dum(79),ilast,iunit
      integer :: idum(4), timeout, istep, npar, srclabel
      integer :: start(500),nwords ! 500 arbitrary
      real    :: time_step, now,xyz(3)
      real    :: amp,xcen,ycen,vdiff,hdiff
      character*256 :: filename
      character*999 :: strbuf, initialization_data

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

c     ---- create all particles once before main time loop ---
c
c     raw particle creation is at a fairly low level because many 
c     variants are anticipated; for simplicity just read one creation site
c
      call read_control_data(simulation_file,"create_particles",strbuf)
      call tokenize(strbuf, start, nwords)
      if      (nwords < 4)  then 
         stop "usage: create_particles = lon,lat,depth,nparticles"
      else if (nwords >= 4) then 
         read(strbuf,*) xyz, npar
         if (nwords == 4) then
            initialization_data = "" 
         else
            initialization_data = strbuf(start(5):) ! rest of line contains bio init data
         endif
      endif
      call setup_ensemble(par_ens, npar) 
      srclabel = 0   ! trace origin of particles
      call generate_particles(par_ens, xyz, time_step, npar, 
     +                        srclabel, initialization_data)
      
      
      current_time => get_master_clock()
c     =====================  main time loop =====================
      istep   = 0
      open(44,file="trajectory1")
      do while (compare_clocks(current_time, end_time)
     +                         *nint(time_step) <= 0)             ! opposite for forward/backward simulation

         write(*,372) istep
         call update_physical_fields()
c        -------- propagate tracers  --------        
         call update_particles(par_ens, time_step)
         
c        -------- write tracer state -------- 
c        -------- loop control       --------
         call add_seconds_to_clock(current_time, nint(time_step))
         istep = istep + 1

         call get_last_particle_number(par_ens, last) 
         if (last>0) then
            call get_particle_position(get_particle(par_ens,1),xyz)
            write(44,*) xyz
         endif
      enddo
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
