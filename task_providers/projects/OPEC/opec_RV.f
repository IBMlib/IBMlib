ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     -----------------------------------------------------------------
c     OPEC phase 1 indices: reproductive volumes (RV and OES conditions)
c     -----------------------------------------------------------------
c     $Rev: $
c     $LastChangedDate:  $
c     $LastChangedBy: $ 
c
c     O2 concentration unit conversion: biological literature often report
c     O2 concentrations as mL/L, where mL is assumed to refer to gas state 
c     (real, not ideal and STP: 1 atm, 0 °C) after extraction from the water. 
c     To convert this to a cencentration measure we apply the USGS recommended conversion 1.42905 mg/mL
c     (ref: Office of Water Quality Technical Memorandum 2011.03,
c     Subject: Change to Solubility Equations for Oxygen in Water)
c     The intrinsic unit in the IMBlib PBI is mmol/m3
c
c     rho calculation code from ~/NSParticleTracking/water_pressure.f
c
c     make ibmrun
c     ibmrun task_providers/projects/OPEC/simpar_RV
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program ibmrun
      use input_parser
      use time_tools
      use run_context
      use physical_fields
      use water_density
      use polygons
      use constants

      implicit none
c     ------------ declarations ------------      
      type(clock),pointer  :: current_time
      type(clock)          :: start_time
      type(clock)          :: end_time
      type(lonlat_polygon), allocatable :: poly(:)
      integer              :: time_step
      real    :: lonmin, lonmax, lon_step
      real    :: latmin, latmax, lat_step
      real    :: ver_step
      character(len=999) :: ofname      
      real    :: xyz(3), cwd, r3(3), dvol
      integer :: idum4(4), ipo, iunit, istat, npoly, ihit
      integer :: nwords, ndigits, nxscan, nyscan, nscan, i, ix,iy
      integer, allocatable :: nstatpt(:), nwetpt(:)
      integer, allocatable :: neu_bouy_faults(:), neu_bouy_wcol(:) 
      integer, allocatable :: neu_bouy_surf(:), neu_bouy_bott(:)
      real,    allocatable :: xyscan(:,:), area(:)
      real,    allocatable :: oesarea(:), reprovol(:)
      real,    allocatable :: totalvol(:), totalarea(:)
      logical, allocatable :: incl_point(:), xymember(:,:)
      
      real           :: nodebuf(1000), BB(4), areafac, lam, phi
      real           :: loctemp, locsalt, locoxy, neu_bouy_z, oes_here
      logical        :: RV_OK
      integer        :: start(1000), ntimpt
      character*999  :: strbuf, pname
      real           :: RV_salt_limit, rho_egg ! read from input
      

      real, parameter :: O2_g_per_L   =  1.42905                       ! USGS, Office of Water Quality Technical Memorandum 2011.03
      real, parameter :: O2_g_per_mol = 15.9994*2                      ! standard atomic weight
      real, parameter :: mll_2_mmolm3 = 1.e3 * O2_g_per_L/O2_g_per_mol ! derived factor for converting ml/l    -> mmol/m3 
      real, parameter :: mmolm3_2_mll = 1.0/mll_2_mmolm3               ! derived factor for converting mmol/m3 -> ml/l
      real, parameter :: RV_oxy_limit  =  2.0*mll_2_mmolm3             ! convert to unit mmol/m3
      real, parameter :: RV_temp_limit =  2.0                          ! unit = degC, ref = Koester unpub
      real, parameter :: neu_bouy_z_resol = 1e-3                       ! unit = m, resolution of neutrally bouyant position

c     ------------   show time starts  ------------
c
      call init_run_context()
c      
c     ------------ set clocks  ------------
c
      call read_control_data(simulation_file, "start_time", idum4)
      call set_clock(start_time, idum4(1), idum4(2), idum4(3), idum4(4))
      write(*,311) "start time = ", idum4
      call read_control_data(simulation_file, "end_time", idum4)
      call set_clock(end_time, idum4(1), idum4(2), idum4(3), idum4(4))
      write(*,311) "end time = ", idum4
      call read_control_data(simulation_file, "time_step", time_step)
      write(*,312) time_step

 311  format(a12,2x,"year",1x,i4,1x,"month",1x,i2,1x,"day",1x,i2,1x,
     +       "second in day",1x,i5,1x)
 312  format("sampling time step = ", 1x,i6,1x,"seconds")
      
      call init_physical_fields(start_time)
      call init_water_density()  ! configure with manual update of plugin
      current_time => get_master_clock()
c
c     ------------ set avaraging polygons + scanning range ------------
c   
      npoly = count_tags(simulation_file, "polygon")  ! expect polygon = name lon1 lat1 lon2 lat2   ...
      allocate( nstatpt(npoly) )
      allocate( nwetpt(npoly)  )
      allocate( poly(npoly)    )
      allocate( oesarea(npoly)   )  ! absolute reproductive area per polygo (by OES condition)
      allocate( reprovol(npoly)  )  ! absolute reproductive volume per polygon (by RV condition)
      allocate( totalvol(npoly)  )  ! total volume per polygon (with same discretization error, to generate volume fractions)
      allocate( totalarea(npoly) )  ! total area per polygon (with same discretization error, to generate area fractions)
      allocate( neu_bouy_faults(npoly) )
      allocate( neu_bouy_wcol(npoly) )
      allocate( neu_bouy_surf(npoly) )
      allocate( neu_bouy_bott(npoly) )

      lonmin =  1.e6
      latmin =  1.e6
      lonmax = -1.e6 
      latmax = -1.e6 
      ihit   = 1     ! handle for scanning polygon definions
      do ipo = 1, npoly
         call read_control_data(simulation_file,"polygon",strbuf,ihit)
         call tokenize(strbuf, start, nwords)
         ndigits = nwords - 1
         if (mod(ndigits,2) == 1) then  !  lon,lat pairs -> even number of ndigits required
            write(*,*) "invalid polygon:", trim(adjustl(strbuf))
            stop
         else   ! start(1) is name
            read(strbuf(start(1):),*) pname
            read(strbuf(start(2):),*) nodebuf(1:ndigits)   ! start(2) indicates start of lon1,lat1,lon2,lat2 ...
         endif
         call init_lonlat_polygon(poly(ipo),
     +                            reshape(nodebuf,(/2, ndigits/2 /)), 
     +                            pname)
         call print_lonlat_polygon(poly(ipo), 6)
         call get_bounding_box(poly(ipo), BB)      ! BB = SWlon, SWlat, NElon, NElat
         lonmin = min(lonmin, BB(1))
         latmin = min(latmin, BB(2))
         lonmax = max(lonmax, BB(3))
         latmax = max(latmax, BB(4)) 
         ihit = ihit + 1   ! progress to next polygon definition
      enddo
      call read_control_data(simulation_file, "longitude_step",
     +                          lon_step)
      call read_control_data(simulation_file, "latitude_step",
     +                          lat_step)  

      write(*,321) "longitude", lonmin, lonmax, lon_step
      write(*,321) "latitude",  latmin, latmax, lat_step

 321  format("Statistics for",a10,"range = ",
     +       f9.3,1x,f9.3,1x,"with spacing",1x,f9.3)
c
c     --- preprocess scanning grid ---
c
      nxscan = 1 + int((lonmax-lonmin)/lon_step)
      nyscan = 1 + int((latmax-latmin)/lat_step)
      nscan  = nxscan*nyscan
      allocate( xyscan(2, nscan)  )
      allocate( area(nscan)       )
      allocate( incl_point(nscan)  ) 
      areafac = (earth_radius*deg2rad)**2 * lon_step * lat_step ! m2
      i = 1
      do iy=1,nyscan
         phi = latmin + (iy-1)*lat_step
         do ix=1,nxscan
            lam = lonmin + (ix-1)*lon_step
            xyscan(1,i)  = lam
            xyscan(2,i)  = phi
            area(i)      = areafac * cos(phi*deg2rad)     ! m2
            incl_point(i) = (horizontal_range_check(xyscan(:,i)).and.
     +                     .not.is_land(xyscan(:,i)))  
            i = i + 1
         enddo
      enddo
      write(*,*) "number scan grid points = ", nscan
      write(*,*) "number of valid points on scan grid = ", 
     +            count(incl_point)
c
      allocate( xymember(nscan, npoly) )
      do ipo = 1, npoly
         call polygon_member_mask(poly(ipo), xyscan(:,:), 
     +                            xymember(:,ipo))   
         write(*,677) ipo, count(xymember(:,ipo))
      enddo
 677  format("scan grid points inside polygon",i3," = ",i5)
c
c     ------------ set vertical ranges ------------   
c 
      call read_control_data(simulation_file, "vertical_step", ver_step)
     
      write(*,332) ver_step
 332  format("averaing vertical step = ",f9.3)
c  
c     ------------ set biological parameters ------------   
c   
      call read_control_data(simulation_file, "RV_salt_limit", 
     +                       RV_salt_limit)
      call read_control_data(simulation_file, "rho_egg", rho_egg)

c  
c     ------------ set dump file ------------   
c      
      call read_control_data(simulation_file,"output_file", ofname)
      call find_free_IO_unit(iunit) 
      open(iunit, file=ofname)
      write(*,*) "writing output to new file:",trim(adjustl(ofname))     
c     
c     ============== main space+time loop ==============   
c
      reprovol  = 0.0    ! absolute reproductive volume per polygon
      oesarea   = 0.0 
      totalvol  = 0.0    ! total volume per polygon
      totalarea = 0.0    ! total area per polygon
      nstatpt   = 0       ! total number of sampling points
      nwetpt    = 0       ! total number of horizontal wet sampling points
      ntimpt    = 0       ! total number of time sampling points
      neu_bouy_faults = 0
      neu_bouy_wcol   = 0
      neu_bouy_surf   = 0
      neu_bouy_bott   = 0

      do while (compare_clocks(current_time, end_time) <= 0)
         call update_physical_fields()
         call update_water_density()   ! module applied in manual mode
         ntimpt = ntimpt + 1   ! all areas
c        ------ horizontal loop ------
         do i = 1, nscan
            if (incl_point(i).and.any(xymember(i,:))) then    ! avoid considering points outside any polygons
               xyz(1:2) = xyscan(:,i)
c              --- check OES condition at this lateral position xyz(1:2) 
c                  currently do not make statistics on neu_bouy_z              
               call evaluate_OES_factor(xyz, neu_bouy_z_resol, rho_egg, 
     +                                  oes_here, neu_bouy_z, istat)   ! 0 < oes_here < 1
c               write(*,*) "NBZ=", neu_bouy_z
               do ipo = 1, npoly
                  if (xymember(i,ipo)) then
                     if (istat < 0) then
                        neu_bouy_faults(ipo) = neu_bouy_faults(ipo)  + 1
                     elseif (istat==0) then 
                        neu_bouy_wcol(ipo) = neu_bouy_wcol(ipo) + 1
                     elseif (istat==1) then 
                        neu_bouy_surf(ipo) = neu_bouy_surf(ipo) + 1
                     elseif (istat==2) then 
                        neu_bouy_bott(ipo) = neu_bouy_bott(ipo) + 1
                     endif      ! istat
                     nwetpt(ipo)    = nwetpt(ipo) + 1
                     totalarea(ipo) = totalarea(ipo) + area(i) 
                     oesarea(ipo)   = oesarea(ipo) + oes_here*area(i) ! oes_here = 0, if istat < 0
                  endif         ! if (xymember
               enddo            ! do ipo 
c              ------ vertical loop (all water column) ------
               xyz(3) = 0.0                ! start at surface
               dvol   = area(i)*ver_step
               call interpolate_wdepth(xyz, cwd, istat)
               do while (xyz(3) <= cwd)  
                  call interpolate_temp(xyz,   loctemp, istat)
                  call interpolate_salty(xyz,  locsalt, istat)
                  call interpolate_oxygen(xyz, locoxy,  istat)
                  RV_OK = OK_for_repro(loctemp,locsalt,locoxy)
c                 ------ second polygon loop ------
                  do ipo = 1, npoly
                     if (xymember(i,ipo)) then
                        nstatpt(ipo)  = nstatpt(ipo) + 1
                        totalvol(ipo) = totalvol(ipo) + dvol
                        if (RV_OK) reprovol(ipo) = reprovol(ipo) + dvol
                     endif      ! if xymember 
                  enddo         ! do ipo       (first polygon loop)
                  xyz(3)  = xyz(3) + ver_step
               enddo            ! do while xyz (vertical loop) 
c              ---                 
            endif               ! if incl_point      
         enddo                  ! do i         (horizontal loop)  
         call add_seconds_to_clock(current_time, time_step)
      enddo                     ! while compare_clocks (time loop)  
       
c
c     ------- report basic statistics to output file -------
c      
      do ipo = 1, npoly
         if (nstatpt(ipo) > 0) then  
c     --- output file ---
            write(iunit,*) ipo,  reprovol(ipo)/ntimpt/1.e9, 
     +                           reprovol(ipo)/totalvol(ipo), 
     +                           oesarea(ipo)/ntimpt/1.e6, 
     +                           oesarea(ipo)/totalarea(ipo)
c        --- stdout ---
            write(*,610)
            call  get_polygon_name(poly(ipo), pname)

            write(*,620) ipo, trim(pname)
            write(*,623) ipo, "reproductive volume = ",
     +                   reprovol(ipo)/ntimpt/1.e9, "km3"
            write(*,623) ipo, "relative reproductive volume = ",
     +                   reprovol(ipo)/totalvol(ipo), ""
            write(*,623) ipo, "reproductive area = ",
     +                   oesarea(ipo)/ntimpt/1.e6, "km2"
            write(*,623) ipo, "relative reproductive area = ",
     +                   oesarea(ipo)/totalarea(ipo), ""
            write(*,641) ipo, "in water column",neu_bouy_wcol(ipo) 
            write(*,641) ipo, "at surface",neu_bouy_surf(ipo) 
            write(*,641) ipo, "at bottom", neu_bouy_bott(ipo) 
            write(*,641) ipo, "undeterminedy", neu_bouy_faults(ipo) 
            write(*,625) "total number of time samplings", ntimpt
            write(*,625) "total number of sampling points", nstatpt(ipo)
            write(*,626) "avarage number of horizontal wet points",
     +                1.0*nwetpt(ipo)/ntimpt + 1e-7
            write(*,626) "average number of vertical sampling points", 
     +                1.0*nstatpt(ipo)/nwetpt(ipo) + 1e-7
         else                   ! statistics undefined, write nothing to output
            write(*,624) ipo
         endif
      enddo ! do ipo
      write(*,610)
c
      close(iunit)
      call close_run_context()
      call close_physical_fields()
      write(*,*) "normal end data extraction"
 610  format(90("-"))
 620  format("polygon",1x,i5,1x,"name =",1x,a)     
 623  format("polygon",1x,i5,1x,a35,1x,e12.6,1x,a)
 624  format("polygon",1x,i5,1x,
     +        ": no wet sampling points in range specifications")     
 625  format(a45,1x,i12)
 626  format(a45,1x,f12.2)     
 641  format("polygon",1x,i5,1x,": samplings bouyancy",1x,a,1x,i7)



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      contains
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      logical function OK_for_repro(temp,salt,oxy)
c     --------------------------------------------------
c     Implement the Koster et al 2005 boleean 
c     local condition for reproductivity:
c 
c         mu(oxy >= RV_oxy_limit  &  salt >= RV_salt_limit & temp >= RV_temp_limit)  
c     --------------------------------------------------
      real, intent(in) :: temp,salt,oxy
      OK_for_repro = ((oxy  >= RV_oxy_limit)  .and.
     +                (salt >= RV_salt_limit) .and.
     +                (temp >= RV_temp_limit)) 
c      OK_for_repro = (oxy >= RV_oxy_limit)
c      OK_for_repro = (salt >= RV_salt_limit)
      end function OK_for_repro



      real function homeostatic_OES_probability(temp,salt,oxy)
c     --------------------------------------------------
c     Koester 2005 (fig. 4) parameterization relative viability of
c     cod eggs in terms of dissolved O2:
c
c          OES = (1 - exp(-0.71*O2)**11.63
c
c     where O2 is dissolved oxygen in unit ml/l.
c     Currently no temp,salt dependence ...
c     --------------------------------------------------
      real, intent(in) :: temp,salt,oxy ! oxy unit mmol/m3
      real             :: O2_mlperl     ! unit ml/l
      real             :: O2scaled           
      O2_mlperl = oxy*mmolm3_2_mll      ! mmol/m3 -> ml/l
      O2scaled  = 0.71*O2_mlperl
c     --- check range of O2scaled to avoid exceptions for exp or **     
      i    f (O2scaled < 1.0e-2) then
         homeostatic_OES_probability = 0.0
      elseif (O2scaled > 60) then 
         homeostatic_OES_probability = 1.0
      else
         homeostatic_OES_probability = (1 - exp(-O2scaled))**11.63
      endif
      end function homeostatic_OES_probability


  
      subroutine evaluate_OES_factor(xy, zacc, rhotarget, oes, 
     +                               zpos, status)
c     --------------------------------------------------
c     Implement the Koster et al 2005 OES relative probability of local
c     egg survival at horizontal position xy:
c          OES(xy) = p(S(xy,d(xy)), T(xy,d(xy)), O(xy,d(xy))) 
c     where p(S,T,O) is the probability in terms of local hydrographic 
c     conditions and d(xy) is depth at neutral buoyancy for cod eggs:
c               rhotarget = rho(xy, d(xy))
c     where rhotarget is the egg density and rho(xy,z) is the water density
c     This equation is solved to within accuracy zacc on water depth 
c
c     Return where in the water column zpos = d(xy) eggs were found to be
c     neutrally bouyant. Also classify the status as
c         status =  0: egg neutrally bouyant within the water column
c         status =  1: egg floats in the surface
c         status =  2: egg sinks to the bottom
c         status = -1: stratification unstable  (return oes = zpos = 0.0)
c
c     time variability: p -> rate ??
c     --------------------------------------------------
      real, intent(in)     :: xy(:)      ! (lon,lat)
      real, intent(in)     :: zacc       ! absolute accuracy on solution zpos (if it exist)
      real, intent(in)     :: rhotarget  ! floater target density
      real, intent(out)    :: oes        ! 0 < oes < 1, if status >= 0, else set oes=0
      real, intent(out)    :: zpos       ! 0 <= zpos <= CWD, if status >= 0, else set zpos=0
      integer, intent(out) :: status     ! see above
c
      real       :: cwd, rhos, rhob, rhom, geo(3)
      real       :: tnb, snb, onb, z0, z1, rho0, rho1
      integer    :: istat, nsteps, ist
c     --------------------------------------------------
      geo(1:2) = xy(1:2) 
      call interpolate_wdepth(geo, cwd, istat) 
      geo(3) = 0
      call interpolate_rhow(geo, rhos, istat) 
      geo(3) = cwd
      call interpolate_rhow(geo, rhob, istat) 
c
c     Now determine (zpos, status) corresponding to neutral bouyancy
c  
c     --- check and flag unstable water column
      if (rhob < rhos) then
         status = -1
         zpos   = 0.0
         oes    = 0.0
         return  ! no further analysis
      endif
c     --- check for surface floating (but stable water column)
      if (rhotarget < rhos) then
         status = 1
         zpos   = 0.0
c     --- check for sinking to bottom (but stable water column)
      elseif (rhotarget > rhob) then
         status = 2
         zpos   = cwd     
c
c     Now it is established that: rhos < rhotarget < rhob
c     Locate rho_egg = rho(xy, zpos) by plain bisection
c     Bracket is (z0,rho0) (z1,rho1) 
      else
         nsteps = 1 + int(log(cwd/zacc)/log(2.0))
         z0 =   0.0
         z1   = cwd
         rho0 = rhos
         rho1 = rhob
         do ist = 1, nsteps
            geo(3) = 0.5*(z0+z1)
            call interpolate_rhow(geo, rhom, istat) 
            if (rhom < rhotarget) then ! advance left bracket
               z0   = geo(3)
               rho0 = rhom
            else                ! decrease right bracket
               z1   = geo(3)
               rho1 = rhom
            endif
         enddo  ! ist
         status = 0
         zpos   = geo(3)
      endif   ! if (rhob ...
      call interpolate_temp(geo,   tnb, istat)
      call interpolate_salty(geo,  snb, istat)
      call interpolate_oxygen(geo, onb,  istat)
      oes = homeostatic_OES_probability(tnb,snb,onb)
      
      end subroutine evaluate_OES_factor

      


      end program
