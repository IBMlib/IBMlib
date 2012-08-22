ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     MvD run
c
c     Unresolved issues: 
c
c     $Rev:  $
c     $LastChangedDate: $
c     $LastChangedBy: $ 
c   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program eulerian_test
      use input_parser
      use time_tools
      use run_context
      use physical_fields
      use advection_diffusion ! Eulerian toolkit
c      use particles         ! Lagrangian
      implicit none

c     ------------ declarations ------------      
      type(clock),target    :: start_time, end_time
      type(clock),pointer   :: current_time
      real                  :: time_step ! Eulerian
      type(eulerian_field)  :: myfield
      real,pointer          :: dat(:,:,:), cbar(:,:)
      real                  :: geoSW(2), geoNE(2)   ! currently no depth is passed
      real                  :: data0(2), dataBC(2) ! buffers
      integer               :: status, idum4(4), istep,navg,nxs,nys
      character(len=12)     :: ofname
      character*(*),parameter :: fmt = '(9999(f8.3,1x))'
c     ------------   show time starts  ------------
      call init_run_context()

c.....set clocks 
      call read_control_data(simulation_file, "start_time", idum4)
      call set_clock(start_time, idum4(1), idum4(2), idum4(3), idum4(4))
      call read_control_data(simulation_file, "end_time", idum4)
      call set_clock(end_time, idum4(1), idum4(2), idum4(3), idum4(4))
      call read_control_data(simulation_file, "time_step", time_step)

      call init_physical_fields(start_time)
      call update_physical_fields()  ! need topology to verify emission boxes


      call init_advection_diffusion()
c
      call read_control_data(simulation_file, "eulerian_box_SW", geoSW)
      call read_control_data(simulation_file, "eulerian_box_NE", geoNE)
      call read_control_data(simulation_file, 
     +                       "eulerian_data_init_value", data0)
      call read_control_data(simulation_file, 
     +                       "eulerian_data_DirichletBC", dataBC)

      call init_eulerian_field(myfield, geoSW, geoNE, data0, dataBC)
c
      current_time => get_master_clock()
      dat => get_data(myfield)  ! index offsets are passed to dat
      nxs = size(dat, 1) - 2 ! interior domain 1:nxs
      nys = size(dat, 2) - 2 ! interior domain 1:nys
      allocate (cbar(0:nxs+1, 0:nys+1)  ) ! include booundary
      navg = 0
      cbar = 0.0

c     =====================  main time loop =====================
      istep   = 0
      open(48,file="consump.dat") ! see eulerian_dynamics/logistic_with_spatial_mortality.f
c
      do while (compare_clocks(current_time, end_time) <= 0)
         write(*,372) istep
c
         call update_physical_fields()
         call update_eulerian_auxillaries()  
c        -------- propagate tracers  --------
         call update_eulerian_field(myfield, time_step, status)
c        -------- write tracer state -------- 
c        -------- loop control       --------
         call add_seconds_to_clock(current_time, nint(time_step))
         istep = istep + 1
c
         dat => get_data(myfield)
         cbar = cbar + dat(:,:,1)
         navg = navg + 1
      enddo
 372  format(25("*"), " main time loop step ",i5.5, " ", 25("*"))

c     ------------- dump final simulation output -------------
      dat(:,:,1) = cbar/navg
      call dump_2Dframe(myfield, 1, "c_avg.xy")
      
c     ----------------- close down ---------------------------
     
      call close_eulerian_field(myfield)
      call close_advection_diffusion()
      call close_physical_fields()
      
      write(*,*) "normal end of simulation"

      

      end program
