ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Compute water flux over a given transect for a given period
c     on a fixed-depth grid (which may include the entire water column)
c     Flux normal is counter-clockwise to the transect vector
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
      use geometry
     
      implicit none
c     ------------ declarations ------------      
      type(clock),target  :: start_time, end_time
      type(clock),pointer :: current_time
      real    :: time_step
      integer :: ih, iz, mh, mz, idum4(4), istep 
      integer :: nstations, inext, ist, istat, iout=89
      character(len=999)  :: filename
      real                :: uvw(3), dlonlat(2), area, r6(6), r2(2)
      real                :: dr(2), drhat(2), hflux, hflux_intg
      real                :: thickness, z0, z1
      real                :: jacob(2), dz, geo(3), dist, tot_area
      real, allocatable   :: xy(:,:), dl(:), nvec(:,:), z(:)

c     ------------   show time starts  ------------
      call init_run_context()
  
c     ------------   set clocks  ------------   
      call read_control_data(simulation_file, "start_time", idum4)
      call set_clock(start_time, idum4(1), idum4(2), idum4(3), idum4(4))
      call read_control_data(simulation_file, "end_time", idum4)
      call set_clock(end_time, idum4(1), idum4(2), idum4(3), idum4(4))
      call read_control_data(simulation_file, "time_step", time_step)   ! seconds

c     ------------   init system  ------------   
      call init_physical_fields(start_time)
      current_time => get_master_clock()
      call update_physical_fields()

c     ------------  read transect grid   ------------  

      call read_control_data(simulation_file,"transect",  r6)  ! lon0[deg]  lat0[deg]  depth0[m] lon1[deg]  lat1[deg]  depth1[m]
      call read_control_data(simulation_file,"resolution",r2)  ! dx[m]  dz[m]  (horizontal, vertical)

c     --- horizontal sampling ---

      call get_horizontal_distance(r6(1:2), r6(4:5), dist) ! width of bounding box in meters
      mh = max(1, int(0.5 + dist / r2(1))) ! positive definite
      allocate( xy(2,mh) )   ! horizontal sampling points (face centers)
      allocate( dl(mh) )     ! length of faces corresponding to sampling points
      allocate( nvec(3,mh) ) ! Cartesian face normals
      dlonlat  = (r6(4:5) - r6(1:2))/mh  ! dlon,dlat
      nvec = 0.0             ! zero vertical component
      do ih = 1, mh
         xy(:,ih) = r6(1:2) + (ih-0.5)*dlonlat ! face mid point
         call get_xyz2cart_jacobian(xy(:,ih), jacob)
         dr = jacob(1:2)*dlonlat
         dl(ih) = sqrt(sum(dr*dr))
         drhat(1)     = -dr(2)  ! rotate counter clockwise
         drhat(2)     =  dr(1)  ! rotate counter clockwise
         nvec(1:2,ih) =  drhat/sqrt(sum(drhat*drhat)) 
      enddo

c     --- vertical sampling ---

      
      thickness = abs(r6(6)-r6(3))  ! vertical thickness of transect
      z0 = min(r6(6), r6(3))        ! upper point
      z1 = max(r6(6), r6(3))        ! lower point
      mz = max(1, int(0.5 + thickness/r2(2))) ! positive definite
      allocate(  z(mz) )      ! vertical sampling points
      dz = thickness / mz     ! vertical step > 0
      do iz = 1, mz
         z(iz) = z0 + (iz-0.5)*dz ! face mid point
      enddo

      tot_area = thickness*dist ! of the full rectangular transect
      write(*,*) "Scanning transect from ", r6(1:2), " to ", r6(4:5)
      write(*,*) "Depth range from ", z0, " to ", 
     +           z1, " m below ses surface"
      write(*,*) "Horizontal sampling     = ", mh
      write(*,*) "Vertical sampling       = ", mz
      write(*,*) "Transect area           = ", tot_area, " m2"
      write(*,*) "Average transect normal = ", sum(nvec, dim=2)/mh

      call read_control_data(simulation_file,"outputfile", filename)  !
      open(iout, file=filename)
      write(iout, 232) "# time[days] after start", 
     +                 "avg transport length [km]", 
     +                 "avg transect current [m/s]" 
 232  format(a30,   2(5x, a30))
 233  format(f30.6, 5x, f30.6, 5x, e30.6)     

c   
c     =====================  main time loop =====================
c
      istep      = 1
      hflux_intg = 0.0
      do while (compare_clocks(current_time, end_time) <= 0)
         call add_seconds_to_clock(current_time, nint(time_step)) ! loop apply to time step end point 
         call update_physical_fields()
         hflux = 0.0 ! zero accumulant
         do iz = 1, mz
            geo(3) = z(iz)
            do ih = 1, mh       ! horizontal mesh loop
               geo(1:2) = xy(:,ih)
               call interpolate_currents(geo, uvw, istat)
               if (istat > 0) uvw = 0.0 ! Faulty interpolations (whatever reason) are padded with 0
               area = dl(ih)*dz ! m2, area element associated with this sampling
               hflux = hflux + sum(nvec(:,ih)*uvw)*area ! unit m3/sec
c               write(*,*) iz, ih, uvw, sum(nvec(:,ih)*uvw)
            enddo
         enddo
c      
         hflux_intg = hflux_intg + hflux*time_step ! integrated water flux, unit m3
c     
         
         write(*,   265) istep,  hflux_intg,  hflux 
         write(iout,233) istep*time_step/86400.0, 
     +                   hflux_intg/tot_area/1.0d3, 
     +                   hflux/tot_area   ! unit: days, km, m/s,
         istep = istep + 1
      enddo
      close(iout)

 265  format("STEP", 1x, i8, 2x, "FLUXINTG = ", e14.6, " m3   FLUX = ", 
     +       e14.6, " m3/s")

      write(*,*) "normal end of simulation"

      end program
