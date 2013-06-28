ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Read a plain list of horizontal space+time sampling points (stations)
c     and generate a corresponding file of 
c
c         station_number   temperature(surf,bott,avg)  salinity(surf,bott,avg)
c
c     station number comes from list of sampling points
c     ---------------------------------------------------
c     $Rev: 147 $
c     $LastChangedDate: 2010-11-19 00:33:57 +0100 (Fri, 19 Nov 2010) $
c     $LastChangedBy: mpay $ 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program ibmrun
      use input_parser
      use time_tools
      use run_context
      use physical_fields
      
     
      implicit none
c     ------------ declarations ------------ 
      type(clock),target  :: time
      integer       :: idum4(4),istat,iver,iland
      integer       :: npt, ipt, statnum
      real          :: geo(3),uvw(3),temp,wdepth,vstep,ss,tt
      real          :: tsurf,tbott,tavg, ssurf,sbott,savg
      logical       :: isla, first
      character*999 :: fname

c     ------------   show time starts  ------------
      call init_run_context()
c
      call read_control_data(simulation_file,"stations", fname)
      open(34,file=fname)
      call read_control_data(simulation_file,"output", fname)
      open(35, file=fname, action='write')
      call read_control_data(simulation_file,"vertical_step", vstep)
c
c     ------ loop over stations ------
c
      ipt      = 0    ! processed points
      iland    = 0    ! land points (discarded)
      first    = .true.
      do
         read(34,*,end=555) statnum, idum4, geo(1:2)
         call set_master_clock(idum4)
         if (first) then
           call init_physical_fields()
           first = .false.
         endif
         call update_physical_fields()
         if (is_land(geo)) then
            iland = iland + 1
            write(*, 824) statnum
            cycle
         else        
            call interpolate_wdepth(geo,wdepth,istat)
            geo(3) = wdepth ! surface
            call interpolate_temp(geo,tbott,istat)
            call interpolate_salty(geo,sbott,istat)
            geo(3) = 0 ! surface
            call interpolate_temp(geo,tsurf,istat)
            call interpolate_salty(geo,ssurf,istat)
c           --- vertical loop ---          
            tavg   = tsurf
            savg   = ssurf
            iver   = 1
            geo(3) = vstep
            do while (geo(3)<wdepth)
               call interpolate_temp(geo,tt,istat)
               call interpolate_salty(geo,ss,istat)
               tavg = tavg + tt
               savg = savg + ss
               iver = iver + 1
               geo(3) = geo(3) + vstep
            enddo
            tavg = tavg / iver  ! iver >= 1
            savg = savg / iver  ! iver >= 1
            write(35,823) statnum, tsurf,tbott,tavg, ssurf,sbott,savg
            ipt = ipt + 1
         endif
      enddo

 823  format(i5, 2x, 3f9.3, 2x, 3f9.3)
 824  format("station", 1x, i5, 1x, "was dry - proceeding to next") 

 555  close(35)
      call close_physical_fields()
      write(*,*) "processed ", ipt, " stations ->", trim(adjustl(fname))
      write(*,*) iland, " dry stations was skipped"
      write(*,*) "normal end of run"
     
      end program

