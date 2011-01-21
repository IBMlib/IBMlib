      program ibmrun
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Test template for physical fields interface
c     
c     Below it will "measure" physical properties in the point xyz
c     specified in the input file, e.g. test_suite/test_pbi/simpar
c     This is a good initial test of a new interface
c   
c     1) Point PHYSICAL_FIELDS_DIR to the new interface in config.mk
c     2) Copy it to IBMlib base directory
c     3) Build & run in IBMlib base directory:   
c          make ibmrun
c          ibmrun test_suite/test_pbi/simpar
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
      real    :: xyz(3),lld(3),res3(3),res2(2),res
      integer :: istat, idum4(4)
      
c     ------------   initialise  ------------
      call init_run_context()
c.....set clocks 
      call read_control_data(simulation_file, "start_time", idum4)
      call set_clock(start_time, idum4(1), idum4(2), idum4(3), idum4(4))
      call set_master_clock(start_time)

c     Init fields      
      call init_physical_fields()
      call update_physical_fields()
      call read_control_data(simulation_file, "xyz", xyz)
      

c     ------------   test interface  ------------
      write(*,*) "-----------------------------"
      write(*,*) "Interpolated values"
      write(*,*) "-----------------------------"
      write(*,*) "Position         : ",xyz
      call interpolate_currents(xyz,res3,istat)
      write(*,*) "Currents         : ",res3
      call interpolate_turbulence(xyz,res3,istat)
      write(*,*) "Turbulence       : ",res3
      call interpolate_turbulence_deriv(xyz,res3,istat)
      write(*,*) "Turbulence deriv : ",res3
      call interpolate_wdepth(xyz,res,istat)
      write(*,*) "Depth            : ",res
      write(*,*) "Is_land          : ",is_land(xyz)
      write(*,*) "Is_wet           : ",is_wet(xyz)
      call interpolate_temp(xyz,res,istat)
      write(*,*) "Temp             : ", res
      call interpolate_salty(xyz,res,istat)
      write(*,*) "Salinity         :", res
c      call interpolate_wind(xyz,res2,istat)
c      write(*,*) "Wind             :", res2
c      write(*,*) "-----------------------------"      
c     ------------   end  ------------      
      write(*,*) "normal end of simulation"
      call close_physical_fields()
      end program
