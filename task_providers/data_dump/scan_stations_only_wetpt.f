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
      integer       :: nstations, nlines
      integer       :: idum4(4),istat,idum3(3),i
      integer       :: npt, ipt, statnum
      real          :: geo(3),uvw(3),wdepth, geo_ask(2),dt
      real          :: tbott, sbott, obott, ubott, nh4bott
      real          :: no3bott, po4bott, pombott, chlbott
      real          :: t_min, s_min, o_min, u_min, nh4_min
      real          :: no3_min, po4_min, pom_min, chl_min
      real          :: t_avg, s_avg, o_avg, u_avg, nh4_avg
      real          :: no3_avg, po4_avg, pom_avg, chl_avg
      real          :: t_max, s_max, o_max, u_max, nh4_max
      real          :: no3_max, po4_max, pom_max, chl_max
      logical       :: first, island
      character*999 :: fname, wetptname
      character*32  :: stat_name
      character*32  :: nstats
      integer       :: year_range(2), isa, last_year,im
      integer       :: nsamp, isamp,iland
      real, allocatable :: wetpt(:,:)
   
c     ------------   show time starts  ------------
      call init_run_context()
c
      call read_control_data(simulation_file,"stations", fname)
      open(1,file=fname)
      nlines = 0
      do
         read(1,*,end=10)
         nlines = nlines + 1
      end do
 10   close(1)

      nstations = nlines - 1
      print*, ""
      print*, "There are ", nstations, " stations"
      
      open(34,file=fname)

      call read_control_data(simulation_file,"output", fname)
      open(35, file=fname, action='write')
      write(35,822) "station",
     +     "isLand", "lon", "lat", "extr_depth", "depth", 
     +     "t_min", "t_avg", "t_max",
     +     "s_min", "s_avg", "s_max",
     +     "o_min", "o_avg", "o_max",
     +     "u_min", "u_avg", "u_max",
     +     "nh4_min", "nh4_avg", "nh4_max",
     +     "no3_min", "no3_avg", "no3_max",
     +     "po4_min", "po4_avg", "po4_max",
     +     "pom_min", "pom_avg", "pom_max",
     +     "chl_min", "chl_avg", "chl_max"       
      call read_control_data(simulation_file,"samplings_per_day", nsamp)
      dt = 86400./(1+nsamp)

      call init_physical_fields()

      call read_control_data(simulation_file,"wet_point_file", 
     +     wetptname)
      wetptname = adjustl(wetptname)
      i = index(wetptname,"*")
      if (i<1) then
         write(*,*) "input @ wet_point_file: name should contain *"
         stop 63
      endif
c     ----- replace "*" with grid selection      
      wetptname = wetptname(1:i-1)//trim(data_set_id)//wetptname(i+1:) ! data_set_id exported from PBI for FishHab
      
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
         island = is_land(geo_ask)
         if (island) then
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
         nh4_min = 1e6
         no3_min = 1e6
         po4_min = 1e6
         pom_min = 1e6
         chl_min = 1e6
         t_avg = 0.0
         s_avg = 0.0
         o_avg = 0.0
         u_avg = 0.0
         nh4_avg = 0.0
         no3_avg = 0.0
         po4_avg = 0.0
         pom_avg = 0.0
         chl_avg = 0.0
         t_max = -1e6
         s_max = -1e6
         o_max = -1e6
         u_max = -1e6
         nh4_max = -1e6
         no3_max = -1e6
         po4_max = -1e6
         pom_max = -1e6
         chl_max = -1e6
c        ----- intra day loop -----
         do isamp = 1, nsamp
            idum4(4) = isamp*dt
            call set_master_clock(idum4)
            call update_physical_fields()
c                   
            call interpolate_wdepth(geo,wdepth,istat)
c     Take the value at the bottom (or very close to it)
            geo(3) = wdepth-0.5 
            call interpolate_temp(geo,tbott,istat)
            call interpolate_salty(geo,sbott,istat)
            call interpolate_oxygen(geo,obott,istat)
            call interpolate_currents(geo,uvw,istat)
            call interpolate_nh4(geo,nh4bott,istat)
            call interpolate_no3(geo,no3bott,istat)
            call interpolate_po4(geo,po4bott,istat)
            call interpolate_part_org_matter(geo,pombott,istat)
            print*, pombott
            call interpolate_chlorophyl(geo,chlbott,istat)
            ubott = sqrt(sum(uvw**2))
            t_avg = t_avg + tbott
            s_avg = s_avg + sbott
            o_avg = o_avg + obott
            u_avg = u_avg + ubott
            nh4_avg = nh4_avg + nh4bott
            no3_avg = no3_avg + no3bott
            po4_avg = po4_avg + po4bott
            pom_avg = pom_avg + pombott
            chl_avg = chl_avg + chlbott
            t_min = min(t_min, tbott)
            s_min = min(s_min, sbott)
            o_min = min(o_min, obott)
            u_min = min(u_min, ubott)
            nh4_min = min(nh4_min, nh4bott)
            no3_min = min(no3_min, no3bott)
            po4_min = min(po4_min, po4bott)
            pom_min = min(pom_min, pombott)
            chl_min = min(chl_min, chlbott)
            t_max = max(t_max, tbott)
            s_max = max(s_max, sbott)
            o_max = max(o_max, obott)
            u_max = max(u_max, ubott)
            nh4_max = max(nh4_max, nh4bott)
            no3_max = max(no3_max, no3bott)
            po4_max = max(po4_max, po4bott)
            pom_max = max(pom_max, pombott)
            chl_max = max(chl_max, chlbott)
         enddo ! isamp
         t_avg = t_avg / nsamp
         s_avg = s_avg / nsamp   
         o_avg = o_avg / nsamp
         u_avg = u_avg / nsamp
         nh4_avg = nh4_avg / nsamp
         no3_avg = no3_avg / nsamp
         po4_avg = po4_avg / nsamp
         pom_avg = pom_avg / nsamp
         chl_avg = chl_avg / nsamp
c        ----- dump data for this station -----
         write(35,823) trim(stat_name),
     +        island, geo(1), geo(2), geo(3), wdepth,
     +        t_min, t_avg, t_max,
     +        s_min, s_avg, s_max,
     +        o_min, o_avg, o_max,
     +        u_min, u_avg, u_max,
     +        nh4_min, nh4_avg, nh4_max,
     +        no3_min, no3_avg, no3_max,
     +        po4_min, po4_avg, po4_max,
     +        pom_min, pom_avg, pom_max,
     +        chl_min, chl_avg, chl_max
c
      enddo    ! station loop

 822  format(a32, 2x, a9, x, 4a12, 9(3a9))
 823  format(a32, 2x, l9, x, 4f12.6, 9(3f9.3))
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

      
