      program ibmrun
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Dump variables for physical fields interface
c     
c     Below it will "dump" one set of physical properties in the surface 
c     layer between the points xy1 and xy2 specified in the input file, 
c     e.g.  test_suite/dump_pbi/2d_field/simpar
c     This is a good initial test of a new interface
c   
c     1) Point PHYSICAL_FIELDS_DIR to the new interface in config.mk
c     2) Copy it to IBMlib base directory
c     3) Build & run in IBMlib base directory:   
c          make ibmrun
c          ibmrun test_suite/dump_pbi/2D_field/simpar
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use input_parser
      use time_tools
      use run_context
      use physical_fields
      use particles
     
      implicit none
c     ------------ declarations ------------      
      type(clock),target  :: start_time
      real    :: xy1(2), xy2(2),xyz(3),dx,dy,res3(3),res,dxy(2) 
      integer :: n_lat,n_lon, i,j,istat,idum4(4),iunit 
      character(99) :: output_var
      
c     ------------   initialise  ------------
      call init_run_context()
c.....set clocks 
      call read_control_data(simulation_file, "start_time", idum4)
      call set_clock(start_time, idum4(1), idum4(2), idum4(3), idum4(4))
      call set_master_clock(start_time)

c     Init fields      
      call init_physical_fields()
      call update_physical_fields()
      call read_control_data(simulation_file, "xy1", xy1)
      call read_control_data(simulation_file, "xy2", xy2)
      call read_control_data(simulation_file, "output_var", output_var)
      call read_control_data(simulation_file, "n_lon", n_lat)
      call read_control_data(simulation_file, "n_lat", n_lon)
      output_var = adjustl(output_var)   ! Get rid of hangers on
      call toupper(output_var)

c     ------------   test interface  ------------
      write(*,*) "============================="
      write(*,*) trim(output_var) //" field"
      write(*,*) "-----------------------------"
      call find_free_IO_unit(iunit)
      open(unit=iunit,file="2D_fields.txt")
      write(iunit,*) "lon lat " // trim(output_var)
      dxy=xy2-xy1
      xyz(3) =0
      do i=1,n_lat
        do j=1,n_lon
          xyz(1) = xy1(1) + dxy(1)*(i-1)/(n_lat-1)
          xyz(2) = xy1(2) + dxy(2)*(j-1)/(n_lon-1)
          select case (output_var(1:4))
          case ("TEMP")
            call interpolate_temp(xyz,res,istat)
          case ("SALI")
            call interpolate_salty(xyz,res,istat)
          case ("U   ")
            call interpolate_currents(xyz,res3,istat)
            res = res3(1)
          case ("V   ")
            call interpolate_currents(xyz,res3,istat)
            res = res3(2)
          case ("HDIF")
            call interpolate_turbulence(xyz,res3,istat)
            res = res3(1)
          case ("VDIF")
            call interpolate_turbulence(xyz,res3,istat)
            res = res3(3)
          case ("DEPT")
            call interpolate_wdepth(xyz,res,istat)
          case default
            call abort_run("dump_pbi_field","Unknown output variable"
     +        // " requested : " //trim(output_var) )
          end select
          write(iunit,*) xyz(1:2), res
        enddo
      enddo
      close(iunit)
      write(*,*) "-----------------------------"
      call close_physical_fields()
      end program
