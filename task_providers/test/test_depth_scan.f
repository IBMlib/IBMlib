ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Test template for physical fields interface
c     ---------------------------------------------------
c     $Rev: 147 $
c     $LastChangedDate: 2010-11-19 00:33:57 +0100 (Fri, 19 Nov 2010) $
c     $LastChangedBy: mpay $ 
c
c     make ibmrun
c     ibmrun task_providers/test/testpar
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program ibmrun
      use input_parser
      use time_tools
      use run_context
      use physical_fields
      use particles
     
      implicit none
c     ------------ declarations ------------      
      integer      :: idum4(4),istat
      integer      :: npt, ipt
      type(clock),target  :: start_time
      integer :: idum(4), timeout, istep, ntracers
      real    :: geo(3),uvw(3),temp,wdepth
      logical :: isla
c     ------------   show time starts  ------------
      call init_run_context()

c.....set clocks 
      call read_control_data(simulation_file, "start_time", idum4)
      call set_clock(start_time, idum4(1), idum4(2), idum4(3), idum4(4))
      call init_physical_fields(start_time)
      call update_physical_fields()
c
      write(*,*) "start_time = "
      call write_clock(start_time)
      call read_control_data(simulation_file,"horizontal_pos", geo(1:2))
      write(*,*) "horizontal_pos          = ", geo(1:2)
      isla = is_land(geo)
      write(*,*) "is_land(horizontal_pos) = ", isla
      if (isla) stop "please select a wet point"
      
      call interpolate_wdepth(geo,wdepth,istat)
      call check_interpol(istat,"wdepth@horizontal_pos")
      write(*,*) "wdepth@horizontal_pos = ", wdepth, "m"

      call read_control_data(simulation_file,"n_interpol_pts", npt)
      write(*,*) "vertical profile with ", npt, "points"
      write(*,423) "depth[m]", "uvw[m/s]", "temp[C]"
      do ipt = 0,npt
         geo(3) = ipt*wdepth/npt
         call interpolate_currents(geo,uvw,istat)
         call check_interpol(istat,"currents")
         call interpolate_temp(geo,temp,istat)
         call check_interpol(istat,"temperature")
         write(*,422) geo(3),uvw, temp
      enddo
 422  format(f8.3, 2x, 3f8.3, 2x, f8.3)
 423  format(  a8, 2x,   a24, 2x,   a8)     
c
      call close_physical_fields()
      
      contains

      subroutine check_interpol(i,msg)
      integer, intent(in)       :: i
      character*(*), intent(in) :: msg
      if (i /= 0) then
         write(*,*) "interpolation",msg,"failed with exit status",i
      endif
      end subroutine check_interpol

      end program

