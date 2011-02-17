
      program bee_test
c     ------------------------------------------------------------------------
c     Compile & run:  ifort -e90 coast_line_intersection_hash_tables_general_mesh.f; a.out
c     ------------------------------------------------------------------------
c     Plot results in R       
c     dat <- read.table("fort.10");
c     plot(NA,asp=1,xlim=range(dat[,c(1,3,5,7)]),ylim=range(dat[,c(2,4,6,8)]));
c     segments(dat$V1,dat$V2,dat$V3,dat$V4);
c     segments(dat$V5,dat$V6,dat$V7,dat$V8);
c     abline(v=-180:180+0.5,h=-90:90+0.5,col="lightgrey",lty=3)
c     ----------------------------------------------------------------
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
      real    :: stability
      type(spatial_attributes) :: bee
      integer :: istep, write_int
c     ----------------------------------------------------------------

c     Initialise modules
      write(*,*) "Bee test:  $Rev: 258 $"
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
      

c     Initialise the bee
      call set_tracer_position(bee,start_geo)
     
      do istep = 1, nint(10**stability)
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

      
