ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     At a specified point advance time
c     and assess currents, temperature and salinity
c
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
      type(clock),target  :: start_time, end_time
      type(clock),pointer :: current_time
      integer :: idum(4), timeout, istep, ntracers
      real    :: geo(3),uvw(3),temp,wdepth,salt,time_step
      real    :: reldepth
      logical :: isla
c     ------------   show time starts  ------------
      call init_run_context()

c.....set clocks 
      call read_control_data(simulation_file, "start_time", idum4)
      call set_clock(start_time, idum4(1), idum4(2), idum4(3), idum4(4))
      call read_control_data(simulation_file, "end_time", idum4)
      call set_clock(end_time, idum4(1), idum4(2), idum4(3), idum4(4))
      call read_control_data(simulation_file, "particle_time_step", 
     +     time_step)
      
      call init_physical_fields(start_time)
      call update_physical_fields()
c
      write(*,*) "start_time = "
      call write_clock(start_time)
      write(*,*) "end_time = "
      call write_clock(end_time)
      
      call read_control_data(simulation_file,"horizontal_pos", geo(1:2))
      write(*,*) "horizontal_pos          = ", geo(1:2)
      isla = is_land(geo)
      write(*,*) "is_land(horizontal_pos) = ", isla
      if (isla) stop "please select a wet point"
 
      call interpolate_wdepth(geo,wdepth,istat)
      call check_interpol(istat,"wdepth@horizontal_pos")
      write(*,*) "wdepth@horizontal_pos = ", wdepth, "m"

      call read_control_data(simulation_file,"relative_depth", reldepth)
      geo(3) = reldepth*wdepth
      write(*,*) "sampling depth@horizontal_pos = ", geo(3), "m"
      
      write(*,423) "time[days]", "uvw[m/s]", "temp[C]","salt[PSU]"

c     =====================  main time loop =====================
      current_time => get_master_clock()
      istep   = 0
      do while (compare_clocks(current_time, end_time) <= 0)
         
         call interpolate_currents(geo,uvw,istat)
         call check_interpol(istat,"currents")
         call interpolate_temp(geo,temp,istat)
         call check_interpol(istat,"temperature")
         call interpolate_salty(geo,salt,istat)
         call check_interpol(istat,"salinity")
         write(*,422) istep*time_step/86400., uvw, temp, salt
         call add_seconds_to_clock(current_time, nint(time_step))
         call update_physical_fields()
         istep = istep + 1
      enddo

c      call dump_land_points("jj.ixiy")
      
 422  format(f10.3, 2x, 3f10.3,      2(2x,f10.3))
 423  format(  a10, 2x, 10x,a10,10x, 2(2x,a10))     
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

