ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     -----------------------------------------------------------------
c     OPEC phase 1 indices
c     -----------------------------------------------------------------
c     $Rev: $
c     $LastChangedDate:  $
c     $LastChangedBy: $ 
c
c     rho calculation code from ~/NSParticleTracking/water_pressure.f
c
c     make ibmrun
c     ibmrun task_providers/data_dump/simpar_seasonal_indices
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program ibmrun
      use input_parser
      use time_tools
      use run_context
      use physical_fields
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
      real    :: depthmin, depthmax, ver_step
      character(len=999) :: ofname      
      real    :: xyz(3), cwd, stop_depth, r3(3)
      integer :: idum4(4), ipo, iunit, istat, npoly, ihit
      integer :: nwords, ndigits, nxscan, nyscan, nscan, i, ix,iy, iz
      integer, allocatable :: nstatpt(:), nwetpt(:), ntimpt(:)
      real,    allocatable :: temp(:), xyscan(:,:), area(:)
      real,    allocatable :: oesarea(:), reprovol(:)
      real,    allocatable :: totalvol(:), totalarea(:)
      logical, allocatable :: incl_point(:), xymember(:,:)
      real           :: nodebuf(1000), BB(4), areafac, lam, phi
      real           :: loctemp, locsalt, locoxy
      real           :: tvec(1000), svec(1000), oxvec(1000)  ! vertical scan buffers -> OES evaluation      
      integer        :: start(1000)
      character*999  :: strbuf, pname

      real, parameter :: RV_oxy_limit  = ???
      real, parameter :: RV_salt_limit = 11.0 ! unit = PSU, ref = Koester 2005
      east/west diff ??
      real, parameter :: rho_egg = ??

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
      current_time => get_master_clock()
c
c     ------------ set avaraging polygons + scanning range ------------
c   
      npoly = count_tags(simulation_file, "polygon")  ! expect polygon = name lon1 lat1 lon2 lat2   ...
      allocate( nstatpt(npoly) )
      allocate( nwetpt(npoly)  )
      allocate( ntimpt(npoly)  ) 
      allocate( poly(npoly)    )
      allocate( temp(npoly)    )    ! temperature accumulant
      allocate( oesarea(npoly)   )  ! absolute reproductive area per polygo (by OES condition)
      allocate( reprovol(npoly)  )  ! absolute reproductive volume per polygon (by RV condition)
      allocate( totalvol(npoly)  )  ! total volume per polygon (with same discretization error, to generate volume fractions)
      allocate( totalarea(npoly) )  ! total area per polygon (with same discretization error, to generate area fractions)
      
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
      call read_control_data(simulation_file, "vertical_range", r3)
      depthmin = r3(1)
      depthmax = r3(2)
      ver_step = r3(3)
      write(*,332) depthmin, depthmax, ver_step
 332  format("averaing vertical range = ",f9.3,1x,f9.3,1x,
     +       "with spacing",1x,f9.3)
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
      temp      = 0.0    ! temperature accumulant
      reprovol  = 0.0    ! absolute reproductive volume per polygon
      oesarea   = 0.0 
      totalvol  = 0.0    ! total volume per polygon
      totalarea = 0.0    ! total area per polygon
      nstatpt   = 0       ! total number of sampling points
      nwetpt    = 0       ! total number of horizontal wet sampling points
      ntimpt    = 0       ! total number of time sampling points

      do while (compare_clocks(current_time, end_time) <= 0)
         call update_physical_fields()
         ntimpt = ntimpt + 1   ! all areas
c        ------ horizontal loop ------
         do i = 1, nscan
            if (incl_point(i).and.any(xymember(i,:))) then    ! avoid considering points outside any polygons
               xyz(1:2) = xyscan(:,i)
               call interpolate_wdepth(xyz, cwd, istat)   
               stop_depth = min(depthmax, cwd)
               
               do ipo = 1, npoly
                  if (xymember(i,ipo)) then
                     nwetpt(ipo) = nwetpt(ipo) + 1
                     totalarea(ipo) = totalarea(ipo) + area(i) 
                  endif
               enddo ! do ipo 
c              ------ vertical loop ------
               xyz(3) = depthmin
               iz     = 0 ! apply prior increment
               do while (xyz(3) <= stop_depth)  
                  call interpolate_temp(xyz,   loctemp, istat)
                  call interpolate_salty(xyz,  locsalt, istat)
                  call interpolate_oxygen(xyz, locoxy,  istat)
                  iz = iz+1 ! prior increment
                  tvec(iz)  = loctemp  ! add to vertical buffer
                  svec(iz)  = locsalt  ! add to vertical buffer
                  oxvec(iz) = locoxy   ! add to vertical buffer
c              ------ first polygon loop ------
                  do ipo = 1, npoly
                     if (xymember(i,ipo)) then
                        temp(ipo)    = temp(ipo) + loctemp 
                        nstatpt(ipo) = nstatpt(ipo) + 1
                        dvol         = area(i)*ver_step
                        totalvol(ipo) = totalvol(ipo) + dvol
                        if (OK_for_repro(loctemp,locsalt,locoxy)) then
                           reprovol(ipo) = reprovol(ipo) + dvol
                        endif
                     endif      ! if xymember 
                  enddo         ! do ipo       (first polygon loop)
                  xyz(3)  = xyz(3) + ver_step
c                 --- pass active buffer segments 1:iz to evaluate_OES_factor
                  oes_here = evaluate_OES_factor(tvec(1:iz), 
     +                                           svec(1:iz), 
     +                                           oxvec(1:iz))
                  do ipo = 1, npoly
                     if (xymember(i,ipo)) then
                        oesarea(ipo)  = oesarea(ipo) +  oes_here*area(i)
                     endif      ! if xymember 
                  enddo         ! do ipo       (second polygon loop) 
               enddo            ! do while xyz (vertical loop) 
            endif               ! if incl_point      
         enddo                  ! do i         (horizontal loop)  
         call add_seconds_to_clock(current_time, time_step)
      enddo                     ! while compare_clocks (time loop)  
       
c
c     ------- report basic statistics to output file -------
c      
      do ipo = 1, npoly
         if (nstatpt(ipo) > 0) then  
            temp(ipo) = temp(ipo) / nstatpt(ipo) 
c     --- output file ---
            write(iunit,*) ipo,  temp(ipo)
c        --- stdout ---
            write(*,610)
            write(*,623) ipo, temp(ipo)
            write(*,625) "total number of sampling points", nstatpt(ipo)
            write(*,626) "avarage number of horizontal wet points",
     +                1.0*nwetpt(ipo)/ntimpt(ipo) + 1e-7
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
 623  format("polygon",1x,i5,1x," : spatio-temporal average =",1x,f12.7)
 624  format("polygon",1x,i5,1x,
     +        ": no wet sampling points in range specifications")     
 625  format(a45,1x,i12)
 626  format(a45,1x,f12.2)     
    

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      contains
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      logical function OK_for_repro(temp,salt,oxy)
c     --------------------------------------------------
c     Implement the Koster et al 2005 boleean 
c     local condition for reproductivity:
c 
c         mu(oxy > RV_oxy_limit & salt > RV_salt_limit)  
c     --------------------------------------------------
      real, intent(in) :: temp,salt,oxy
      OK_for_repro = ((oxy >= RV_oxy_limit).and.(salt >= RV_salt_limit))
      end function OK_for_repro


      real function homeostatic_OES_probability(temp,salt,oxy)
c     --------------------------------------------------
c     Wieland 19XX parameterization  
c     --------------------------------------------------
      real, intent(in) :: temp,salt,oxy
      homeostatic_OES_probability = ...
      end function homeostatic_OES_probability



      real function evaluate_OES_factor(tvec, svec, ovec)
c     --------------------------------------------------
c     Implement the Koster et al 2005 OES probability of local
c     egg survival at horizontal position x:
c          OES(x) = p(S(x,d(x)), T(x,d(x)), O(x,d(x))) 
c     where p(S,T,O) is the probability in terms of local hydrographic 
c     conditions and d(x) is depth at neutral buoyancy for cod eggs:
c               rho_egg = rho(x, d(x))
c     where rho_egg = ??? and rho(x,z) is the water density
c
c     time variability: p -> rate ??
c     --------------------------------------------------
      real, intent(in) :: tvec(:), svec(:), ovec(:)
      real             :: rho( size(tvec) ) ! automatic array
      integer          :: iv, nv
      nv = size(tvec) ! assume equal for svec and oxvec

      SETUP RHO
      SOLVE RHO=RHO_EGG -> Z
      tnb,snb,onb =   
      evaluate_OES_factor = homeostatic_OES_probability(tnb,snb,onb)
      end function evaluate_OES_factor

      


      end program
