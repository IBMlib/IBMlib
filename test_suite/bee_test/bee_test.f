ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     -----------------------------------------------------------
c     Bee_test 
c     -----------------------------------------------------------
c     $Rev$
c     $LastChangedDate$
c     $LastChangedBy$ 
c
c     Provides a rigorous, brute-force test of the coastline
c     intersection routine contained in a full working PBI. The 
c     test is predicated on a particle walking randomly around the
c     wet domain space and interacting with the coastline. Options
c     are configured by the external controlfile
c    
c     Original version by ASC. Adapted to work as a task_provider
c     by MPA 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program bee_test
      use run_context
      use input_parser
      use physical_fields
      use particle_tracking
      use particles
      implicit none
      real    :: start_geo(3), deg2rad
      real    :: dgeo(3),geo(3)
      logical :: surfenc, botenc, shoreenc, domainenc
      logical :: write_traj 
      real    :: angle,r2(2)
      real    :: rjump   ! Bee jump range
      real    :: stability,frac
      type(spatial_attributes) :: bee
      integer :: write_int
      integer(kind=8) :: disp_freq,loops,istep
c     ----------------------------------------------------------------

c     Initialise modules
      write(*,*) "Bee test:  $Rev$"
      call init_run_context()
      call init_physical_fields()
      call init_particles()

c     Read from the ctrlfile
      call read_control_data(simulation_file, "stability", stability)
      call read_control_data(simulation_file, "start_geo", start_geo)
      call read_control_data(simulation_file, "step_size", rjump)
      call read_control_data(simulation_file, "write_traj", write_int)
      write_traj = .FALSE.
      if(write_int.ne.0) then
        write(*,*) "Writing trajectory...."
        write_traj=.TRUE.
      endif
      disp_freq = nint(10**(stability-3),8)  
      loops = nint(10**stability,8) 

c     Initialise the bee
      call set_tracer_position(bee,start_geo)
     
      do istep = 1, loops 
         if(mod(istep,disp_freq)==0) then
            frac= real(istep)/real(loops)*100
            write(*,'(a,f5.1,a)') "Bee test ",frac,"% complete."
         endif
         !Setup the step
         call random_number(r2)
         angle = 8*atan(1.0)*r2(2)
         dgeo  = r2(1)*rjump*(/cos(angle), sin(angle), 0.0/)

         !Make the step
         call add_constrained_step(bee,dgeo,
     +            surfenc, botenc, shoreenc, domainenc)
         call get_tracer_position(bee,geo)

         !Write to file if desired
         if(write_traj) then
            write(65,*) geo(1:2)
         end if

c        !check the step
         if (is_land(geo)) then
            write(*,*) "geo dry at step",istep
            write(*,*) "geo = ", geo
            write(*,*) geo(1:2)
           stop
         endif
         
      enddo

      write(*,*) "=========normal end of simulation========="


      end program

      
