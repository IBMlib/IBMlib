ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Minimal template for generating a hydrographic data dump 
c     
c     This example scans a regular longitude-latitude grid
c     and dumps horizontal current in the middle of the 
c     water column
c     ---------------------------------------------------
c     $Rev: $
c     $LastChangedDate:  $
c     $LastChangedBy: $ 
c
c     make ibmrun
c     ibmrun task_providers/data_dump/simpar
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program ibmrun
      use input_parser
      use time_tools
      use run_context
      use physical_fields
     
      implicit none
c     ------------ declarations ------------      
      
      type(clock),target  :: start_time
      integer :: ix,iy,iz,istat
      real    :: xyz(3)
      real    :: lon,lat, uvw(3),wd
      character*256 :: filename
      integer, parameter :: nlon = 20
      integer, parameter :: nlat = 20
      real, parameter    :: lonmin = -8
      real, parameter    :: lonmax = 8
      real, parameter    :: latmin = 52
      real, parameter    :: latmax = 58 
c     ------------   show time starts  ------------
      call init_run_context()
      
c.....set clocks 
      
      call set_clock(start_time, 2004, 1, 5, 0)
      call init_physical_fields(start_time)
      call update_physical_fields()

      xyz = 0.
      open(77,file="uvdat")
      do iy = 1,nlat
         lat = latmin + iy*(latmax-latmin)/nlat
         do ix = 1,nlon
            lon = lonmin + iy*(lonmax-lonmin)/nlon
            xyz(1:2) = (/lon,lat/)
            if (.not.is_land(xyz)) then
               call interpolate_wdepth(xyz,wd,istat)
               xyz(3) = 0.5*wd
               call interpolate_currents(xyz,uvw,istat)
               write(77,*) lon, lat, uvw(1:2)
            else
               write(77,*) lon, lat, 0, 0
            endif
         enddo
      enddo
      write(*,*) "normal end of simulation"

      end program
