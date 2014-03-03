ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     scan tests of si interpolation
c     ---------------------------------------------------
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program ibmrun
      use input_parser
      use time_tools
      use run_context
      use physical_fields
      use atlantis_grid
c     ------------ declarations ------------      
      implicit none
      type(clock),target  :: start_time, end_time
      real                :: geo(3), cwd, dz, dx,value
      integer             :: i, istat,idum4(4)
      integer,parameter   :: nscan = 1000
c     ----------------------------------------------------------------------------------------------------

c     ------------   show time starts  ------------
      call init_run_context()
c  
c.....set clocks 
c      
      call read_control_data(simulation_file,"start_time", idum4)
      call set_clock(start_time, idum4(1), idum4(2), idum4(3), idum4(4))
     
      call init_physical_fields(start_time)
      call update_physical_fields()
c
c.....vertical scan
c
      geo(1:2) = (/3.0, 54.0/)
      call interpolate_wdepth(geo, cwd, istat)
      dz = cwd/nscan
      do i=1,nscan
         geo(3) = i*dz
         call interpolate_salty(geo, value, istat)
         if (istat == 0) then
            write(44,*) geo(3),value
         else
            write(*,*) "vertical scan ", geo(3), istat
         endif
      enddo
c
c.....West-East scan
c
      geo = (/-2.0, 55.0, 0.05/)
      dx = 14.0/nscan
      do i=1,nscan
         geo(1) = -2 + i*dx
         call interpolate_salty(geo, value, istat)
         if (istat /= 0) value = 0.0
         write(45,*) geo(1), value
      enddo
c            
      write(*,*) "normal end of simulation"
c      
      end program
