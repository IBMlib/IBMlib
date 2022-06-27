ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ------------------------------------------------------------------------------------------------------
c     Experimental DRRS (dynamical renormalization resampling scheme) particle distribution equilibrator
c
c     Based on CLAIM macro-plastic work simulating an open source-sink system
c     ------------------------------------------------------------------------------------------------------
c     This version only supports wet-boundary condition c = 0 and is based on projection on a finite source set
c
c     Next step will be c=const, but this requires defining an effective advective velocity, since also
c     stokes + wind drift contribute. Further, vertical effects needs to be adressed consistently.
c     Final option could be dynamic wet-boundary condition c = clocal(t) or c = avgerage(t)
c
c
c     This task module requires an extension of the particle state API, namely the function with signature:
c          real function get_sink_rate(state_attribute)
c
c     Further the task module requires an extension of physics API (only needed by write_aux_data, can be removed):
c          subroutine interpolate_wind_10m(xyz, uv, istat)
c          subroutine interpolate_surf_stokes(xyz, uv, istat)
c
c     ------------------------------------------------------------------------------------------------------
c     $Rev$
c     $LastChangedDate$
c     $LastChangedBy$ 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program tracker
      use write_netcdf_lonlat_data  
      use input_parser
      use time_tools
      use run_context
      use physical_fields
      use particle_state
      use particles
      use particle_tracking
      use geometry
      use array_tools
     
      implicit none
c     ------------ declarations ------------      
      integer      :: idum4(4)
      type(clock),target  :: start_time, end_time
      type(clock),pointer :: current_time
      type(particle_ensemble)     :: par_ens
      type(emission_box), pointer :: emitboxes(:)
      type(netcdf_grid_2Ddata_handler) :: nch            ! for output
      type(netcdf_grid_2Ddata_handler) :: nch1,nch2,nch3 ! for optional output
      integer :: year, month, day, julday,last,nxo,nyo,iofreq
      integer :: i,ix,iy,iz,dum(79),ilast,iunit,inew,ic,ihit,nobl
      integer :: idum(4), timeout, istep, ntracers, nreshuf,istat
      real    :: time_step, now,xyz(3), srcseed, reshuf_rate
      real    :: amp,xcen,ycen,vdiff,hdiff, risk, tmem, alpha
      real    :: outputfreq, pos(3),wd
      character*256 :: filename
      real, allocatable :: wet_src_pts(:,:)   ! allowed (lon,lat) for sources; shape = nwpts,2
      real, allocatable :: influx(:)          ! per sources; shape = nwpt
      real, allocatable :: src_chances(:)     ! chance of receiving an agent; shape = nwpt
      real, allocatable :: my_luck(:)         ! lucky straws for particles; shape = npart
      real, allocatable :: open_boundary_line(:,:) ! shape (ic,4)
      real, allocatable :: area_km2(:,:)      ! output grid cells
      real, parameter   :: sec_in_year = 86400.0*365.0
      logical           :: write_aux
      type(spatial_attributes), pointer :: all_spaces(:)
      type(state_attributes), pointer   :: all_states(:)
c     ------------   show time starts  ------------
      call init_run_context()

c.....set clocks 
      call read_control_data(simulation_file, "start_time", idum4)
      call set_clock(start_time, idum4(1), idum4(2), idum4(3), idum4(4))
      call read_control_data(simulation_file, "end_time", idum4)
      call set_clock(end_time, idum4(1), idum4(2), idum4(3), idum4(4))
      call read_control_data(simulation_file, "particle_time_step", 
     +                       time_step)
      call init_physical_fields(start_time)
      call init_particles()
      call update_physical_fields() ! need topology to verify emission boxes
      
c     ---- topography available now  ----
      
      call setup_ensemble_from_file(par_ens, "emitbox", emitboxes)

      call read_control_data(simulation_file,"wet_point_file", filename)
      call load_wetpt(filename)

      call read_control_data(simulation_file,"source_point_file",
     +                       filename)
      call load_source_influx(filename)

      call read_control_data(simulation_file,"reshuf_memory_time", tmem) ! tmem in days
      alpha = exp(-time_step/tmem/86400.) ! for reshuffling running average
      
      allocate( my_luck( get_ensemble_size(par_ens) ))
c
c     parse open boundary configuration
c
      nobl = count_tags(simulation_file, "open_boundary_line")
      if (nobl > 0) then
         allocate( open_boundary_line(nobl,4) )
         ihit = 1
         do ic = 1, nobl
            call read_control_data(simulation_file,"open_boundary_line",
     +                             open_boundary_line(ic,:),ihit)
            write(*,593) ic, open_boundary_line(ic,:)
         enddo
      endif
 593  format("open_boundary_line ",i2," = ", 4f12.5)

      call init_nc_griddata_handler(nch)   ! calls init_netcdf_dataset
      call init_nc_aux_handler(nch1,nch2,nch3,write_aux)   ! calls init_netcdf_dataset
      call read_control_data(simulation_file,"output_frequency", ! unit days
     +                       outputfreq)
      iofreq = int(0.5 + outputfreq*86400/time_step)
                       
c     =====================  main time loop =====================
      
      current_time => get_master_clock()
      istep   = 0
      
      do while (compare_clocks(current_time, end_time)
     +                         *nint(time_step) <= 0)             ! opposite for forward/backward simulation
         write(*,372) istep
         call update_physical_fields()
c        -------- propagate tracers  --------        
         call generate_particles(par_ens, emitboxes, time_step)
         call update_particles(par_ens, time_step)
c        -------- reshuffle --------
         call random_number(my_luck)
         nreshuf = 0            ! this time step
         all_states => get_active_particle_states(par_ens)
         all_spaces => get_active_spatial_parts(par_ens)
         if (associated(all_states)) then
            last = size(all_states)
         else
            last = 0
         endif
         
         do i=1, last
            risk = get_sink_rate(all_states(i))*time_step
            call get_tracer_position(all_spaces(i), xyz)
            if (violate_confinement(xyz)) risk = 1.0
c           --- tmp rooting out of threspassers (21 Oct 2018)
            if ((.not.horizontal_range_check(xyz))
     +           .or.is_land(xyz)
     +           .or.(.not.is_wet(xyz))) risk = 1.0
c           --- reshuffle unluckly ones             
            if (my_luck(i) < risk) then
               call init_state_attributes(all_states(i), all_spaces(i),
     +                                    time_step, "", 0) ! delete reference to emit box
c              ----- select new source for this particle
               call random_number(srcseed)
               call search_sorted_list(srcseed, src_chances, inew)
               xyz(1:2) = wet_src_pts(inew,1:2) ! transfer source position
               call set_tracer_position(all_spaces(i), xyz)
               nreshuf = nreshuf + 1
            endif
         enddo
         
         if (istep == 0) reshuf_rate = nreshuf                 ! initial condition
         reshuf_rate = alpha*reshuf_rate + nreshuf*(1.0-alpha) ! per time step
         
         write(*,373) istep, nreshuf, last, reshuf_rate

         if (mod(istep, iofreq)==1) then
            call write_distrib(nch, reshuf_rate*sec_in_year/time_step)
            if (write_aux) call write_aux_data(nch1, nch2, nch3)
         endif
c        -------- loop control       --------
         call add_seconds_to_clock(current_time, nint(time_step))
         istep = istep + 1

      enddo

 372  format(25("*"), " main time loop step ",i5.5, " ", 25("*"))
 373  format("step", 1x, i5.5, 1x,"reshuffled",1x,i10,"/",i10,1x,
     +       "particles, running average = ", f12.2)
 377  format(f6.2, 1x, 3f8.3, 1x, i2, 1x, f8.3,i3)


c     ----------------- close down ---------------------------

      call close_netcdf_dataset(nch)
      if (write_aux) call close_aux(nch1, nch2, nch3)
      call close_physical_fields()
      call close_particles()

      write(*,*) "normal termination of simulation"

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      contains
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      

      
      subroutine load_wetpt(fname)
c     -------------------------------------------------------
c     Load wet point allowed as sources from fname
c     File format:
c        lon1 lat1
c        lon2 lat2
c     Validate points are wet - topography must be initialized
c     stop if dry/outside domain points are detected
c     -------------------------------------------------------    
      integer :: ipt, nwpts
      character*(*) :: fname
      real          :: geo(2)
c     -------------------------------------------------------    
      open(75,file=fname)
c     ---- first coumt number of points in file ----
      nwpts = 0
      do
         read(75,*,end=534)
         nwpts = nwpts+1
      enddo
 534  rewind(75)
c     ---- read wet points and validate them ----     
      allocate( wet_src_pts(nwpts,2) )
      do ipt=1,nwpts
         read(75,*) geo
         if (.not.horizontal_range_check(geo)) then
            write(*, 573) ipt, geo
            stop
         endif
         if (is_land(geo)) then
            write(*, 574) ipt, geo
            stop
         endif
         wet_src_pts(ipt,:) = geo ! inside domain and wet
      enddo
      
      close(75)
      write(*,569) nwpts, trim(adjustl(fname))
      
 569  format("load_wetpt: loaded ", i, " wet point list from", a)
 573  format("load_wetpt: point", i, " (",2f8.3,") is outside domain")
 574  format("load_wetpt: point", i, " (",2f8.3,") is dry")
      
      end subroutine load_wetpt

      

      subroutine load_source_influx(fname)
c     -----------------------------------------------------------------
c     load sources :                  influx
c     set influx dstribution vector : src_chances
c     -----------------------------------------------------------------
      character*(*) :: fname
      real    :: apos(2), flux
      integer :: inear,ipt,nwpts, nsrc
      real    :: totflux
c     -----------------------------------------------------------------
      open(75,file=fname)
      nwpts = size(wet_src_pts,1)
      allocate( influx(nwpts) )
      influx = 0.0
      nsrc   = 0
      read(75,*,end=431) ! skip header
      do 
         read(75,*,end=431) apos, flux
         call select_nearest_wetpt(apos, inear)
         influx(inear) = influx(inear) + flux
         nsrc = nsrc + 1
      enddo
 431  close(75)
      if (nsrc<1) then
         write(*,565) trim(adjustl(fname))
         stop
      else
         write(*,566) nsrc, trim(adjustl(fname))
      endif
 565  format("load_source_influx: no sources in", a)     
 566  format("load_source_influx: loaded ", i8," sources from ", a)
      
      
c     ----- set influx dstribution vector : src_chances
      totflux = sum(influx)
      allocate( src_chances(nwpts) )
      src_chances(1) = 0.0
      do ipt = 2, nwpts
         src_chances(ipt) = src_chances(ipt-1) + influx(ipt-1)/totflux
      enddo
c      do ipt = 1, nwpts
c         write(31,*) src_chances(ipt) = src_chances(ipt-1) + influx(ipt-1)/totflux
c         write(32,*) ipt, src_chances(ipt)
c         write(32,*) ipt, influx(ipt)
c      enddo
      end subroutine load_source_influx

      

      subroutine select_nearest_wetpt(apos, inear)
c     -------------------------------------------------------
c     select nearest wet point from a list wetpt(npts, 2)
c     -------------------------------------------------------    
      real, intent(in)     :: apos(:)
      integer, intent(out) :: inear
      integer :: ipt
      real    :: dist, d
c     -------------------------------------------------------
      dist  = 1.0e12 ! far away
      inear = -1
      do ipt = 1, size(wet_src_pts, 1)
         call get_local_distance(apos, wet_src_pts(ipt,:), d)
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
      end subroutine select_nearest_wetpt


      logical function violate_confinement(pos)
c     ---------------------------------------------------------
c     check whether any horizontal open boundary lines are violated
c     by point pos
c     open_boundary_line = (x0,y0, vx0,vy0)
c     (x0,y0) = basept  (vx0,vy0) = normal toward forbidden zone
c     ---------------------------------------------------------
      real, intent(in) :: pos(:)
      real             :: dotp, dxy(2)
      integer          :: ic
c     ----------------------------------------------      
      violate_confinement = .false.
      if (allocated(open_boundary_line)) then
         do ic = 1, size(open_boundary_line,1)
            dxy  = pos(1:2) - open_boundary_line(ic, 1:2)
            dotp = sum(dxy * open_boundary_line(ic, 3:4))
            if (dotp>0) then ! we only need sign
               violate_confinement = .true.
               exit
            endif
         enddo
      endif
      end function violate_confinement


      subroutine init_nc_griddata_handler(nch)
c     -------------------------------------------------------
c     -------------------------------------------------------
      type(netcdf_grid_2Ddata_handler), intent(out) :: nch
      character*256        :: outfilename
      real                 :: xy1(2), dxy(2)
      integer              :: idum2(2),  year, month, day
      character*64         :: toffs_attr
c     -------------------------------------------------------
      call read_control_data(simulation_file,"outputfile", outfilename)
      call read_control_data(simulation_file,"output_grid_dims", idum2)
      call read_control_data(simulation_file,"output_grid_SWcorner",xy1)
      call read_control_data(simulation_file,"output_grid_spacings",dxy)
      call get_time_components(start_time, year, month, day)
      call make_time_attribute(toffs_attr, 'days', year, month, day)
      call init_netcdf_dataset(nch, outfilename,  
     +     idum2(1), idum2(2), xy1(1), xy1(2),
     +     dxy(1), dxy(2), 'xyt', toffs_attr)
c      
      allocate( area_km2(idum2(1), idum2(2)) )
      call get_area_matrix(nch, area_km2)
c
      end subroutine init_nc_griddata_handler


      
      subroutine init_nc_aux_handler(nch1,nch2,nch3,write_aux)
c     -------------------------------------------------------
c     -------------------------------------------------------
      type(netcdf_grid_2Ddata_handler), intent(out) :: nch1,nch2,nch3
      logical, intent(out)                          ::  write_aux
      integer :: m
      character*256        :: fstem
      real                 :: xy1(2), dxy(2)
      integer              :: idum2(2),  year, month, day
      character*64         :: toffs_attr
c     -------------------------------------------------------
      m = count_tags(simulation_file, "write_aux")
      write_aux = .false.
      if (m>0) then
         call read_control_data(simulation_file,"write_aux", fstem)
         write_aux = .true.
      else
         return
      endif
c     --- use same grid
      call read_control_data(simulation_file,"output_grid_dims", idum2)
      call read_control_data(simulation_file,"output_grid_SWcorner",xy1)
      call read_control_data(simulation_file,"output_grid_spacings",dxy)
      call get_time_components(start_time, year, month, day)
      call make_time_attribute(toffs_attr, 'days', year, month, day)
      call init_netcdf_dataset(nch1, trim(fstem)//"_curr.nc",  
     +     idum2(1), idum2(2), xy1(1), xy1(2),
     +     dxy(1), dxy(2), 'xyt', toffs_attr)
      call init_netcdf_dataset(nch2, trim(fstem)//"_wind.nc",  
     +     idum2(1), idum2(2), xy1(1), xy1(2),
     +     dxy(1), dxy(2), 'xyt', toffs_attr)
      call init_netcdf_dataset(nch3, trim(fstem)//"_stok.nc",  
     +     idum2(1), idum2(2), xy1(1), xy1(2),
     +     dxy(1), dxy(2), 'xyt', toffs_attr)
      end subroutine init_nc_aux_handler
      
      subroutine close_nc_griddata_handler(nch)
c     -------------------------------------------------------
c     -------------------------------------------------------
      type(netcdf_grid_2Ddata_handler), intent(inout) :: nch
      call close_netcdf_dataset(nch)
      end subroutine close_nc_griddata_handler

      
      subroutine close_aux(nch1, nch2, nch3)
c     -------------------------------------------------------
c     -------------------------------------------------------      
      type(netcdf_grid_2Ddata_handler), intent(inout) :: nch1,nch2,nch3
      call close_netcdf_dataset(nch1)
      call close_netcdf_dataset(nch2)
      call close_netcdf_dataset(nch3)
      end subroutine close_aux


      
      subroutine write_distrib(nch, reshuf_per_year)
c     -------------------------------------------------------
c     Cast particle distribution into a density field
c     on regular grid (lon, lat) save field as netCDF
c     using data set handler nch
c     
c     particles outside grid is ignored
c     -------------------------------------------------------
      type(netcdf_grid_2Ddata_handler), intent(inout) :: nch
      real, intent(in)     :: reshuf_per_year ! emsemble avg reshufflings per year
c      
      integer              :: ipt, ix, iy, npt, days_since_start, last
      logical              :: inside
      real                 :: content_per_particle, xy(2),sxy(2), hours
      real                 :: xyz(3)
      integer, allocatable :: particles_in_cells(:,:) ! particle count
      real, allocatable    :: cell_content(:,:) ! flux-renormalized cell content
      type(spatial_attributes), pointer :: all_spaces(:)
c     -------------------------------------------------------
      all_spaces => get_active_spatial_parts(par_ens)
      if (associated(all_spaces)) then
         last = size(all_spaces)
      else
         last = 0
      endif
      
      if (last<1) then
         write(*,*) "save_distribution: no particles, no output"
         return
      endif
     
         
      allocate( particles_in_cells(nch%nxs, nch%nys) )
      allocate( cell_content(nch%nxs, nch%nys)       )
      particles_in_cells = 0
      cell_content       = 0.
      
      do ipt = 1, last
         call get_tracer_position(all_spaces(ipt), xyz)   
         call get_subgrid_indices(nch, xyz(1:2), ix, iy, inside)
c           .... ignore particles outside grid ....         
         if (inside) then
            particles_in_cells(ix,iy) = particles_in_cells(ix,iy) + 1
         endif
      enddo
      npt = sum(particles_in_cells)
      if (last>0) write(*,781) 100.0*npt/last
      
c     .... conversion factor from particle to content
      content_per_particle = sum(influx)/max(0.5, reshuf_per_year) ! avoid zero division
      cell_content = content_per_particle*particles_in_cells/area_km2  
      call get_period_length_hour(start_time, current_time, hours)
      days_since_start = nint(hours/24.)
c     .... dump data ....
      call write_netcdf_data(nch,cell_content,timeval=days_since_start)
      
      write(*,782) nch%nxs, nch%nys, trim(adjustl(nch%fname))
      
 781  format("save_distribution: ", f8.4,
     +       " % of particles on analysis grid")
 782  format("save_distribution: saved ", i5," x ",i5,
     +     " data points to file ", a)
      
      deallocate( particles_in_cells )
      deallocate( cell_content       )
c
      end subroutine write_distrib

      
      subroutine write_aux_data(nch1, nch2, nch3)
c     -------------------------------------------------------
c     -------------------------------------------------------
      type(netcdf_grid_2Ddata_handler), intent(inout) :: nch1,nch2,nch3
      real, allocatable    :: uc(:,:), uw(:,:), us(:,:)
      real :: xyz(3), uvw(3),auc,auw,aus, hours
      integer :: ist1,ist2,ist3,ix,iy,days_since_start
c     -------------------------------------------------------
      allocate( uc(nch1%nxs, nch1%nys) )
      allocate( uw(nch2%nxs, nch2%nys) )
      allocate( us(nch3%nxs, nch3%nys) )
      uc = -1e20
      uw = -1e20
      us = -1e20
      do ix=1,nch1%nxs
         do iy=1,nch1%nys
            xyz(1) = nch1%lon1 + (ix-1)*nch1%dlon
            xyz(2) = nch1%lat1 + (iy-1)*nch1%dlat
            xyz(3) = 1e-8
            if (is_land(xyz)) cycle
            call interpolate_currents(xyz, uvw, ist1)
            auc = sqrt(sum(uvw*uvw))
            call interpolate_wind_10m(xyz, uvw(1:2), ist2)
            auw = sqrt(sum(uvw(1:2)*uvw(1:2)))
            call interpolate_surf_stokes(xyz, uvw(1:2), ist3)  
            aus = sqrt(sum(uvw(1:2)*uvw(1:2)))
            if (max(ist1,ist2,ist3)>0) cycle  ! accept/reject all together
            uc(ix,iy) = auc
            uw(ix,iy) = auw 
            us(ix,iy) = aus
         enddo
      enddo
      call get_period_length_hour(start_time, current_time, hours)
      days_since_start = nint(hours/24.)
      call write_netcdf_data(nch1,uc,timeval=days_since_start)
      call write_netcdf_data(nch2,uw,timeval=days_since_start)
      call write_netcdf_data(nch3,us,timeval=days_since_start)
      deallocate( uc )
      deallocate( uw )
      deallocate( us )
      end subroutine write_aux_data
      
      end program
