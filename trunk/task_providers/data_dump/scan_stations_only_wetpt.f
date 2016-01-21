ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     The physical model coastline is often slightly different from
c     the real coastline, most often due to finite resolution.
c     This can be a problem for samplings near the coast line.
c     A pragmatic solution to this is to pick the neareast wet point
c     from a pre-defined data set (which could be cell centers for the underlying 
c     physical model grid) rather than a fragile wet point search 
c     algorithm.
c     The example below shows an example how to do this: We want to assess 
c     the hydrography at a set of sampling stations. A sampling station list
c     is loaded, and a pre-defined data set wet_point_file is loaded.
c     If a sampling station is on land (dry), the nearest point from
c     wet_point_file is picked. See subroutine load_wetpt for format of wet_point_file
c
c     Example reading tags:
c       stations          = ../station_list.dat
c       output            = hydrography_at_station_list.dat               
c       samplings_per_day = 24
c       wet_point_file    = wetpoints.lonlat
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program ibmrun
      use input_parser
      use time_tools
      use run_context
      use geometry
      use physical_fields
      
     
      implicit none
c     ------------ declarations ------------ 
      type(clock),target  :: time
      integer       :: idum4(4),istat,idum3(3)
      integer       :: npt, ipt, statnum
      real          :: geo(3),uvw(3),wdepth, geo_ask(2),dt
      real          :: tbott, sbott, obott, ubott
      real          :: t_min, s_min, o_min, u_min
      real          :: t_avg, s_avg, o_avg, u_avg
      real          :: t_max, s_max, o_max, u_max       
      logical       :: first
      character*999 :: fname, wetptname
      character*32  :: stat_name
      integer       :: year_range(2), isa, last_year,im
      integer       :: nsamp, isamp,iland
      real, allocatable :: wetpt(:,:)
     
c     ------------   show time starts  ------------
      call init_run_context()
c
      call read_control_data(simulation_file,"stations", fname)
      open(34,file=fname)

      call read_control_data(simulation_file,"output", fname)
      open(35, file=fname, action='write')
      write(35,822) "station",
     +              "t_min", "t_avg", "t_max",
     +              "s_min", "s_avg", "s_max",
     +              "o_min", "o_avg", "o_max",
     +              "u_min", "u_avg", "u_max"
        
      call read_control_data(simulation_file,"samplings_per_day", nsamp)
      dt = 86400./(1+nsamp)

      call init_physical_fields()

      call read_control_data(simulation_file,"wet_point_file", 
     +                       wetptname)
      call load_wetpt(wetptname)

c
c     ------ loop over stations ------
c     
c     format:
c         station   year   month   day  lon   lat
      
      statnum = 0
      iland   = 0
      read(34,*) ! skip header line
      do
         read(34,*,end=555) stat_name, idum3, geo_ask
         statnum = statnum + 1 ! process station
         idum4(1:3) = idum3
         write(*,*) "station  = ", stat_name
         write(*,*) "date     = ", idum3
         write(*,*) "position = ", geo_ask
         if (is_land(geo_ask)) then
            call select_nearest_wetpt(geo_ask, geo)
            iland = iland + 1
            write(*,*) "position recast to ", geo
         else
            geo(1:2) = geo_ask
         endif
         t_min = 1e6
         s_min = 1e6
         o_min = 1e6
         u_min = 1e6
         t_avg = 0.0
         s_avg = 0.0
         o_avg = 0.0
         u_avg = 0.0
         t_max = -1e6
         s_max = -1e6
         o_max = -1e6
         u_max = -1e6
c        ----- intra day loop -----
         do isamp = 1, nsamp
            idum4(4) = isamp*dt
            call set_master_clock(idum4)
            call update_physical_fields()
c                   
            call interpolate_wdepth(geo,wdepth,istat)
            geo(3) = wdepth-0.5 
            call interpolate_temp(geo,tbott,istat)
            call interpolate_salty(geo,sbott,istat)
            call interpolate_oxygen(geo,obott,istat)
            call interpolate_currents(geo,uvw,istat)
            ubott = sqrt(sum(uvw**2))
            t_avg = t_avg + tbott
            s_avg = s_avg + sbott
            o_avg = o_avg + obott
            u_avg = u_avg + ubott
            t_min = min(t_min, tbott)
            s_min = min(s_min, sbott)
            o_min = min(o_min, obott)
            u_min = min(u_min, ubott)
            t_max = max(t_max, tbott)
            s_max = max(s_max, sbott)
            o_max = max(o_max, obott)
            u_max = max(u_max, ubott)
         enddo ! isamp
         t_avg = t_avg / nsamp
         s_avg = s_avg / nsamp   
         o_avg = o_avg / nsamp
         u_avg = u_avg / nsamp  
c        ----- dump data for this station -----
         write(35,823) stat_name, 
     +                 t_min, t_avg, t_max,
     +                 s_min, s_avg, s_max,
     +                 o_min, o_avg, o_max,
     +                 u_min, u_avg, u_max   
         
c
      enddo    ! station loop

 822  format(a32, 2x, 4(3a9))
 823  format(a32, 2x, 4(3f9.3))
c     ---- no more data ---      
 555  close(35)
      call close_physical_fields()
      write(*,*) "processed ", statnum, " stations ->",
     +           trim(adjustl(fname))
      write(*,*) iland, " dry stations were recast"
      write(*,*) "normal end of run"
    
 
      contains


      subroutine load_wetpt(fname)
c     -------------------------------------------------------
c     Load wet point list from fname
c     File format:
c        number_of_wet_points_to_follow
c        lon1 lat1
c        lon2 lat2
c        ...
c     -------------------------------------------------------    
      integer :: ipt, nwpts
      character*(*) :: fname
c     -------------------------------------------------------    
      open(75,file=fname)
      read(75,*) nwpts
      allocate( wetpt(nwpts,2) )
      do ipt=1,nwpts
         read(75,*) wetpt(ipt,:)
      enddo
      close(75)
      write(*,*) "load_wetpt: loaded wet point list from ",
     +           trim(adjustl(fname))
      end subroutine load_wetpt
      


      subroutine select_nearest_wetpt(apos, take_this)
c     -------------------------------------------------------
c     select nearest wet point from a list wetpt(npts, 2)
c     -------------------------------------------------------    
      real, intent(in)  :: apos(:)
      real, intent(out) :: take_this(:)
      integer :: ipt, inear
      real    :: dist, d
c     -------------------------------------------------------
      dist  = 1.0e12 ! far away
      inear = -1
      do ipt = 1, size(wetpt, 1)
         call get_local_distance(apos, wetpt(ipt,:), d)
         if (d < dist) then ! update best candidate
            dist  = d
            inear = ipt
         endif   
      enddo
      if (inear<0) then
         write(*,*) "select_nearest_wetpt: unexpected exception"
         stop 323
      endif
c
      take_this(1:2) = wetpt(inear,:)
c     
c      write(*,*) "check:target =",apos
c      write(*,*) "check:nearest=",take_this
c      write(*,*) "check:dist   =",dist
c
      end subroutine select_nearest_wetpt

      end program

      
