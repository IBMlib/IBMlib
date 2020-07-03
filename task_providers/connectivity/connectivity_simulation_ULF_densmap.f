cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     -------------------------------------------------------------------------
c     Standard particle tracking simulation 
c     writing output in Unified Lagrangian Format (ULF)
c     as stored in netcdf :
c
c       section 0) premises (input, version, configurations, ...) - not implemented yet
c       section 1) optional topography (for visualization)
c       section 2) optional particle tracks (number of frames bind to unlimited dimension)
c       section 3) optional connectivity matrix                 
c
c     This version implement sections 1+2+3
c     -------------------------------------------------------------------------
c     https://www.unidata.ucar.edu/software/netcdf/netcdf-4/newdocs/netcdf-f90/
c     -------------------------------------------------------------------------
c     add density maps feature to connectivity_simulation_ULF.f 
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program tracker
      use netcdf
      use input_parser
      use time_tools
      use run_context
      use physical_fields
      use particles
      use particle_tracking
      use particle_state     !  get_settlement_habitats
      use connectivity
      use write_netcdf_lonlat_data ! for creating density map
      
      implicit none
c     ------------ declarations ------------      
      type(clock),target  :: start_time, end_time
      type(clock),pointer :: current_time
      type(particle_ensemble)     :: par_ens
      type(emission_box), pointer :: emitboxes(:)                
       
      type(state_attributes), pointer :: state_stack(:)
      integer :: year, month, day, julday
      integer :: i,ix,iy,iz,dum(79),ilast,iunit, numsec, npart,np
      integer :: idum0(4), idum1(4), timeout, istep, ntracers, ntim
      integer :: ncid, nskip, topolon_id, topolat_id, npart_id
      integer :: lon_id, lat_id, depth_id, time_id, i4_id, nframes_id
      integer :: lon0_id, lat0_id, depth0_id
      integer :: nxtop, nytop, nxtop_id, nytop_id, topography_id
      integer :: is_out_of_domain_id, is_settled_id
      integer :: age_at_settling_id, degree_days_id, max_temp_id
      integer :: min_temp_id, max_salt_id, min_salt_id
      integer :: idum2(2)
      real    :: time_step, t, rdum6(6), dt_save
      real, allocatable :: topography(:,:), topolon(:), topolat(:)
      logical :: save_tracks, save_xy, save_xyz, topo_scan
      logical :: save_tracer_stat
      character*256 :: filename, topofilename, sconn, aword
      integer :: save_connectivity
      integer :: nsrc, ndest, nsrc_id, ndest_id, conmat_id
      integer :: nemax
      real, allocatable    :: tmat(:,:), tmat_prior(:)
      integer, allocatable :: num_emit(:)

      integer :: nxs, nys, noutside, i2(2), i4(4), last
      real    :: lon1,lat1, dlon, dlat, pos(3), r4(4), r2(2)
      logical :: inside, incl_part
      INTEGER*8, allocatable :: particle_count(:,:) ! range -9223372036854775808, 9223372036854775807
      character(len=64) :: vattr(2,2) 
      type(netcdf_grid_2Ddata_handler) :: nchout  ! type from write_netcdf_lonlat_data
      real :: min_temp_allowed, max_temp_allowed
      real :: min_salt_allowed, max_salt_allowed
      
c     ------ for windows debugging bind stdout/stderr to specific files (non portable) ------

      open(0, file="stderr.txt")           ! comment to allow shell redirection of stderr
c      open(6, file="stdout.txt")          ! comment to allow shell redirection of stdout

c     ------------   show time starts  ------------
      call init_run_context()

c      
c.....set clocks 
c      
      call read_control_data(simulation_file, "start_time", idum0)
      call set_clock(start_time, idum0(1), idum0(2), idum0(3), idum0(4))
      call read_control_data(simulation_file, "end_time", idum1)
      call set_clock(end_time, idum1(1), idum1(2), idum1(3), idum1(4))
      call read_control_data(simulation_file, "particle_time_step", 
     +                       time_step)
      call init_physical_fields(start_time)
      
      call get_period_length_sec(start_time, end_time, numsec)
      ntim = 1 + int(numsec/time_step) ! number of time steps in this simulation
c
c.....setup ensemble generator      
c
      call init_particles()
      call update_physical_fields() ! need topology to verify emission boxes


      call create_emission_boxes("emitbox",  emitboxes, np)   ! only activated in all time steps
 
      call setup_ensemble(par_ens, np)
      npart = get_ensemble_size(par_ens) ! at this point we know max number of particles
      if (npart<1) then
         write(*,*) "WARNING: you have no particles in your simulation"
      endif
c
c.....check whether vertical hard limits apply
c
      if (count_tags(simulation_file, "hard_vertical_limits") > 0) then
         call read_control_data(simulation_file,
     +        "hard_vertical_limits", idum2)
         par_ens%min_depth = idum2(1)
         par_ens%max_depth = idum2(2)
         write(*,*) "imposing hard_vertical_limits =", idum2
      else
         write(*,*) "no hard_vertical_limits imposed"
      endif
      
c.....prepare output   
c 
      call read_control_data(simulation_file, "outputfile", filename)
      call nfcheck( nf90_create(filename, NF90_CLOBBER, ncid) )
      
c
c.....prepare output: section 1 (topography)    
c     
c     set control handle topo_scan, possibly allocate/define arrays
c
      if (count_tags(simulation_file, "topography_file") > 0) then
         call read_control_data(simulation_file, "topography_file", 
     +                    topofilename)      
         call nfcheck( nf90_put_att(ncid, NF90_GLOBAL, 
     +                 "topography_file", topofilename) )
         topo_scan = .false.
      elseif (count_tags(simulation_file, "sample_topography") > 0) then
         call nfcheck( nf90_put_att(ncid, NF90_GLOBAL, 
     +                 "topography_scan", "yes") )
         topo_scan = .true.
         call read_control_data(simulation_file, "sample_topography",
     +                          rdum6) 
c        rdum6 = lon0, lon1, dlon ; lat0, lat1, dlat
c                  1    2     3       4    5     6
         nxtop = 1 + int( (rdum6(2)-rdum6(1)) / rdum6(3) )
         nytop = 1 + int( (rdum6(5)-rdum6(4)) / rdum6(6) ) 
         allocate ( topography(nytop, nxtop) )
         allocate ( topolon(nxtop)        )
         do ix = 1, nxtop
            topolon(ix)  = rdum6(1) + (ix-1)*rdum6(3)
         enddo
         allocate ( topolat(nytop)         )
         do iy = 1, nytop
            topolat(iy)  = rdum6(4) + (iy-1)*rdum6(6)
         enddo
         call nfcheck( nf90_def_dim(ncid, "nxtop", nxtop, nxtop_id) )
         call nfcheck( nf90_def_dim(ncid, "nytop", nytop, nytop_id) )
c        -------- define topolon --------
         call nfcheck( nf90_def_var(ncid, "topolon", NF90_FLOAT, 
     +                 (/nxtop_id/), topolon_id) )
         call nfcheck( nf90_put_att(ncid, topolon_id, "meaning", 
     +                 "topography sampling zonal grid")) 
         call nfcheck( nf90_put_att(ncid, topolon_id, "unit", 
     +                 "degrees East")) 
c        -------- define topolat --------
         call nfcheck( nf90_def_var(ncid, "topolat", NF90_FLOAT, 
     +                 (/nytop_id/), topolat_id) )
         call nfcheck( nf90_put_att(ncid, topolat_id, "meaning", 
     +                 "topography sampling meridional grid")) 
         call nfcheck( nf90_put_att(ncid, topolat_id, "unit", 
     +                 "degrees North")) 
c        -------- define topography --------
         call nfcheck( nf90_def_var(ncid, "topography", NF90_FLOAT, 
     +                 (/nytop_id, nxtop_id/), topography_id) )
         call nfcheck( nf90_put_att(ncid, topography_id, "meaning", 
     +                 "topography sampling on grid")) 
         call nfcheck( nf90_put_att(ncid, topography_id, "unit", 
     +                              "meters below sea surface"))
c
      else             ! section 1 empty
         topo_scan = .false.
      endif
c
c.....prepare output: section 2 (tracks)    
c     
c     set control handles save_tracks, save_xy, save_xyz, possibly allocate/define arrays
c
      if (count_tags(simulation_file, "save_xyz_frequency") > 0) then
         call read_control_data(simulation_file, "save_xyz_frequency", 
     +                     dt_save)
         save_xy  = .false.
         save_xyz = .true.
      elseif (count_tags(simulation_file, "save_xy_frequency") > 0) then
         call read_control_data(simulation_file, "save_xy_frequency", 
     +                     dt_save)
         save_xy  = .true.
         save_xyz = .false.
      else
         save_xy  = .false.
         save_xyz = .false.
         dt_save  = 1.0  ! but render defined to avoid exeptions
      endif
      save_tracks = save_xy .or. save_xyz
      if (dt_save > 0.0) then
         nskip = max(1, nint(dt_save/time_step)) ! skip frames between saves
      else
         nskip = ntim                            ! save only last time frame
      endif

      save_tracer_stat = .false. ! default
      if (count_tags(simulation_file, "save_tracer_stat") > 0) then
         call read_control_data(simulation_file, "save_tracer_stat",
     +        aword)
         if (trim(adjustl(aword)) == "yes") save_tracer_stat = .true.
      endif
      
c     NF90_UNLIMITED binds to particle frames

      if (save_tracks) then
         call nfcheck( nf90_put_att(ncid, NF90_GLOBAL, 
     +        "tracks_xy", "yes") )
         call nfcheck( nf90_def_dim(ncid, "nframes", NF90_UNLIMITED, 
     +        nframes_id))
         call nfcheck( nf90_def_dim(ncid, "nparticles",
     +        max(1,npart), npart_id))   ! max(1,npart): have at least one slot to avoid clash with NF90_UNLIMITED
         call nfcheck( nf90_def_dim(ncid, "i4",        4,     i4_id))
c        -------- define lon --------
         call nfcheck( nf90_def_var(ncid, "lon", NF90_FLOAT, 
     +                 (/npart_id, nframes_id/), lon_id))
         call nfcheck( nf90_put_att(ncid, lon_id, "meaning", 
     +                 "particles zonal position")) 
         call nfcheck( nf90_put_att(ncid, lon_id, "unit", 
     +        "degrees East"))
         
         call nfcheck( nf90_def_var(ncid, "lon0", NF90_FLOAT, 
     +                 (/npart_id/), lon0_id))
         call nfcheck( nf90_put_att(ncid, lon0_id, "meaning", 
     +                 "particles zonal initial position")) 
         call nfcheck( nf90_put_att(ncid, lon0_id, "unit", 
     +        "degrees East"))
         
c        -------- define lat --------
         call nfcheck( nf90_def_var(ncid, "lat", NF90_FLOAT, 
     +                 (/npart_id, nframes_id/), lat_id))
         call nfcheck( nf90_put_att(ncid, lat_id, "meaning", 
     +                 "particles meridional position")) 
         call nfcheck( nf90_put_att(ncid, lat_id, "unit", 
     +        "degrees North"))
         
         call nfcheck( nf90_def_var(ncid, "lat0", NF90_FLOAT, 
     +                 (/npart_id/), lat0_id))
         call nfcheck( nf90_put_att(ncid, lat0_id, "meaning", 
     +                 "particles meridional initial position")) 
         call nfcheck( nf90_put_att(ncid, lat0_id, "unit", 
     +        "degrees North"))
         
c        -------- define time --------
         call nfcheck( nf90_def_var(ncid, "time", NF90_INT, 
     +                 (/i4_id, nframes_id/), time_id))
         call nfcheck( nf90_put_att(ncid, time_id, "meaning", 
     +                 "time corresponding to particles frames")) 
         call nfcheck( nf90_put_att(ncid, time_id, "unit", 
     +                            "year, month, day, second_in_day"))
      endif
c     ----------------
      if (save_xyz) then
         call nfcheck( nf90_put_att(ncid, NF90_GLOBAL, 
     +                 "tracks_xyz", "yes") )
c        -------- define depth --------
         call nfcheck( nf90_def_var(ncid, "depth", NF90_FLOAT, 
     +        (/npart_id, nframes_id/), depth_id))   
         call nfcheck( nf90_put_att(ncid, depth_id, "meaning", 
     +                 "particles vertical position")) 
         call nfcheck( nf90_put_att(ncid, depth_id, "unit", 
     +        "meters below sea surface"))
         call nfcheck( nf90_def_var(ncid, "depth0", NF90_FLOAT, 
     +        (/npart_id/), depth0_id))   
         call nfcheck( nf90_put_att(ncid, depth0_id, "meaning", 
     +                 "particles vertical position")) 
         call nfcheck( nf90_put_att(ncid, depth0_id, "unit", 
     +        "meters below sea surface"))
         
      endif
c     ---------- save track statistics / final state ----------      
      if (save_tracer_stat) then
         call nfcheck( nf90_def_var(ncid,"is_out_of_domain",NF90_BYTE, 
     +        (/npart_id/), is_out_of_domain_id))
         call nfcheck( nf90_def_var(ncid,"is_settled",      NF90_INT, 
     +        (/npart_id/), is_settled_id))
         call nfcheck( nf90_def_var(ncid,"age_at_settling", NF90_FLOAT,
     +        (/npart_id/), age_at_settling_id))        
         call nfcheck( nf90_def_var(ncid,"degree_days",     NF90_FLOAT,
     +        (/npart_id/), degree_days_id))        
         call nfcheck( nf90_def_var(ncid,"max_temp",        NF90_FLOAT,
     +        (/npart_id/), max_temp_id))        
         call nfcheck( nf90_def_var(ncid,"min_temp",        NF90_FLOAT,
     +        (/npart_id/), min_temp_id))         
         call nfcheck( nf90_def_var(ncid,"max_salt",        NF90_FLOAT,
     +        (/npart_id/), max_salt_id))         
         call nfcheck( nf90_def_var(ncid,"min_salt",        NF90_FLOAT,
     +        (/npart_id/), min_salt_id))         
      endif
      
      
c
c.....prepare output: section 3 (connectivity)    
c     
c     map input tag save_connectivity to integer variable save_connectivity:
c     save_connectivity  =  0: do not save connectivity matrix
c     save_connectivity  =  1: save connectivity matrix corresponding to last time step
c
      save_connectivity = 0 ! also matches absent tag save_connectivity 
      if (count_tags(simulation_file, "save_connectivity") > 0) then
         call read_control_data(simulation_file, "save_connectivity", 
     +        sconn)
         if (trim(adjustl(sconn)) == "last") save_connectivity = 1
      endif

!     --- prepare writing last connectivity matrix ---
      if (save_connectivity == 1) then
         nsrc  = size(emitboxes)  ! emitboxes has been initialized
         ndest = size(get_settlement_habitats()) ! particular particle_state method
         call nfcheck( nf90_put_att(ncid, NF90_GLOBAL, 
     +        "connectivity_final", "yes"))
         call nfcheck( nf90_def_dim(ncid, "nsource", nsrc, nsrc_id))
         call nfcheck( nf90_def_dim(ncid, "ndest", ndest, ndest_id))
         call nfcheck( nf90_def_var(ncid, "connectivity_matrix", 
     +                 NF90_FLOAT, (/nsrc_id, ndest_id/), conmat_id) )
         call nfcheck( nf90_put_att(ncid, conmat_id, "meaning", 
     +                 "probability of particle transport")) 
         call nfcheck( nf90_put_att(ncid, conmat_id, "unit", 
     +                 "0 < probability < 1")) 
      endif

c
c     ==== leave define mode - after here we can write data variables ====
c
      call nfcheck( nf90_enddef(ncid) )  
c      

      if (topo_scan) then
         call scan_topography() ! topolon,topolat -> topography
         call nfcheck( nf90_put_var(ncid, topolon_id, topolon) )
         call nfcheck( nf90_put_var(ncid, topolat_id, topolat) )
         call nfcheck( nf90_put_var(ncid, topography_id, topography) )
      endif

c
c     ------- read density output options -------
c      
      call read_control_data(simulation_file, "density_file", filename)
      call read_control_data(simulation_file, "density_grid_dims",i2)    ! nxs,nys
      nxs = i2(1)
      nys = i2(2)
      call read_control_data(simulation_file, "density_grid_bbox",r4)    ! W/S/E/N 
      dlon = (r4(3)-r4(1))/i2(1)
      dlat = (r4(4)-r4(2))/i2(2)
      lon1 = r4(1) + 0.5*dlon    ! cell center of first grid point
      lat1 = r4(2) + 0.5*dlat   ! cell center of first grid point
      allocate ( particle_count(nxs, nys) )
      particle_count = 0
      vattr(1,1) = "long_name"
      vattr(2,1) = "accumulated particle count"
      vattr(1,2) = "unit"
      vattr(2,2) = "1"
      call init_netcdf_dataset_write(nchout, filename,  
     +     nxs, nys, lon1, lat1, dlon, dlat,
     +     "xy", vname="particle_count",vattr=vattr)
      call read_control_data(simulation_file,
     +     "allowed_temperature_range", r2)
      min_temp_allowed = r2(1)
      max_temp_allowed = r2(2)
      write(*,*) "allowed_temperature_range = ", r2
      
      call read_control_data(simulation_file,
     +     "allowed_salinity_range", r2)
      min_salt_allowed = r2(1)
      max_salt_allowed = r2(2)
      write(*,*) "allowed_salinity_range = ", r2
c      
c     =====================  main time loop =====================
c
      current_time => get_master_clock()
      t = 0.0
      
      do istep = 1, ntim
         write(*,372) istep
         call update_physical_fields()      
         call generate_particles(par_ens, emitboxes,
     +                           time_step) ! activated all time steps
         call update_particles(par_ens, time_step)
         t  = t + time_step
         call add_seconds_to_clock(current_time, nint(time_step))
         if ( save_tracks .and. (mod(istep, nskip) == 0)) then
            call write_particle_frame()   ! corresponding to current_time
         endif
         call write_progress_file(istep,ntim) ! text file for progress monitoring
c        ---- update density map ----
         call get_last_particle_number(par_ens, last)
         do i=1,last
            call get_particle_position(get_particle(par_ens,i), pos)
            call get_subgrid_indices(nchout, pos, ix, iy, inside)
            incl_part = .true.
            if (par_ens%space_stack(i)%outofdomain) incl_part=.false.   ! FTH exclude condition 1a (at boundary)
            if (par_ens%state_stack(i)%type /= particle_free)           ! FTH exclude condition 1b (unfree)
     +          incl_part=.false.   
            if (par_ens%state_stack(i)%bio%min_temp < min_temp_allowed) ! FTH exclude condition 2 
     +           incl_part=.false.
            if (par_ens%state_stack(i)%bio%max_temp > max_temp_allowed) ! FTH exclude condition 2 
     +           incl_part=.false.
            if (par_ens%state_stack(i)%bio%min_salt < min_salt_allowed) ! FTH exclude condition 2 
     +           incl_part=.false.
            if (par_ens%state_stack(i)%bio%max_salt > max_salt_allowed) ! FTH exclude condition 2 
     +           incl_part=.false.
            
c            write(*,*) par_ens%space_stack(i)%outofdomain,
c     +           par_ens%state_stack(i)%type,
c     +           par_ens%state_stack(i)%bio%min_temp,
c     +           par_ens%state_stack(i)%bio%max_temp,
c     +           par_ens%state_stack(i)%bio%min_salt,
c     +           par_ens%state_stack(i)%bio%max_salt
            
            if (inside.and.incl_part) then
               particle_count(ix, iy) = particle_count(ix, iy) + 1
            endif 
         enddo ! i=1,last
      enddo    ! istep = 1, ntim
 372  format(25("*"), " main time loop step ",i5.5, " ", 25("*"))

c     --------- dump density output
      call write_netcdf_data(nchout, real(particle_count, kind=4))  ! API supports real*4
      call close_netcdf_dataset(nchout)
      
c     --------- generate and save final connectivity, if requested ---------

      if (save_connectivity == 1) then
         allocate( tmat(nsrc, ndest) )
         allocate( num_emit(nsrc)    )
         allocate( tmat_prior(nsrc)  )
         do i = 1, nsrc
            call get_current_emissions(emitboxes(i), num_emit(i), nemax)
            tmat_prior(i) = 1.0/ndest/max(1, num_emit(i))
         enddo
         state_stack => get_active_particle_states(par_ens)
         call get_connectivity_matrix(state_stack, num_emit,
     +        tmat_prior, tmat)
         call nfcheck( nf90_put_var(ncid, conmat_id, tmat) )
         write(*,*) "wrote connectivity matrix for last time step"
      endif
c     --------- save track statistics / final state , if requested ---------

      if (save_tracer_stat) call save_tracer_stat_to_netcdf()
      if ( save_tracks )    call write_particle_initial_pos()
      
      
c     ----------------- close down ---------------------------
     
      call nfcheck( nf90_close(ncid) )
      call close_physical_fields()
      call close_particles()
      if (allocated ( topography )) deallocate ( topography )
      if (allocated ( topolon    )) deallocate ( topolon    )
      if (allocated ( topolon    )) deallocate ( topolat    ) 

      write(*,*) "normal termination of simulation"


c     ------------ internal subroutines below  ------------
            contains 
c     ----------------------------------------------------- 

      subroutine nfcheck(status)
      integer, intent ( in) :: status
    
      if(status /= nf90_noerr) then 
         write(0,*)  trim(nf90_strerror(status))
         stop  "Stopped"
      end if

      end subroutine nfcheck        
   

      subroutine scan_topography()
c     -----------------------------------------------------
c     topolon,topolat -> topography(iy,ix)    
c     -----------------------------------------------------
      integer         :: ix,iy,istat
      real            :: pos(2)
      real, parameter :: dry     =  0.0
      real, parameter :: outside = -1.0
c     -----------------------------------------------------      
      do ix = 1, nxtop
         pos(1) = topolon(ix)
         do iy = 1, nytop
            pos(2) = topolat(iy)

            if (horizontal_range_check(pos)) then

               if (is_land(pos)) then
                  topography(iy,ix) = dry
               else
                  call interpolate_wdepth(pos, topography(iy,ix), istat)
               endif

            else

               topography(iy,ix) = outside

            endif

         enddo ! iy
      enddo    ! ix
      end subroutine scan_topography


      subroutine write_progress_file(inow, nsteps)
c     ----- requested by FTH for progress monitoring -----
      integer,intent(in) :: inow, nsteps
      integer            :: iof
      call find_free_IO_unit(iof)
      open(iof, file="SIMSTAT.txt")
      write(iof, 285) "completed steps",     inow
      write(iof, 285) "steps in simulation", nsteps
 285  format(a20, "=",i10)
      close(iof)
      end subroutine write_progress_file

      
      subroutine write_particle_frame()
c     -----------------------------------------------------  
c     reference variables in global scope
c     ----------------------------------------------------- 
      integer             :: i, last, now(4), ifr, dimid
      real                :: xyz(3, npart)
      real, parameter     :: posfill = -1.0
      character*256       :: name
c     ----------------------------------------------------- 
      xyz = posfill
      call get_last_particle_number(par_ens, last) 
      do i=1,last
         call get_particle_position(get_particle(par_ens,i),xyz(:,i))
      enddo
      
      call get_date_from_clock(current_time, now(1), now(2), now(3))
      call get_second_in_day(current_time, now(4))
      call nfcheck( nf90_inquire_dimension(ncid,nframes_id,name,ifr))
      call nfcheck( nf90_put_var(ncid, time_id, now,     
     +              start=(/1, 1+ifr/)) )
      call nfcheck( nf90_put_var(ncid, lon_id,  xyz(1,:),
     +              start=(/1, 1+ifr/)) ) 
      call nfcheck( nf90_put_var(ncid, lat_id,  xyz(2,:),
     +              start=(/1, 1+ifr/)) )
      if (save_xyz) then
         call nfcheck( nf90_put_var(ncid, depth_id,  
     +                 xyz(3,:), start=(/1, 1+ifr/)) )
      endif
      end subroutine write_particle_frame



      
      subroutine write_particle_initial_pos()
c     -----------------------------------------------------  
c     dump initial positions
c     ----------------------------------------------------- 
      integer             :: i, last, now(4), ifr, dimid
      real                :: xyz(3, npart)
      real, parameter     :: posfill = -1.0
      character*256       :: name
c     ----------------------------------------------------- 
      xyz = posfill
      call get_last_particle_number(par_ens, last) 
      do i=1,last
         xyz(:,i) = par_ens%space_stack(i)%init_position
      enddo
      call nfcheck( nf90_put_var(ncid, lon0_id,  xyz(1,:)))
      call nfcheck( nf90_put_var(ncid, lat0_id,  xyz(2,:)))
      if (save_xyz) then
         call nfcheck( nf90_put_var(ncid, depth0_id, xyz(3,:)))
      endif
      end subroutine write_particle_initial_pos
      

      
      subroutine save_tracer_stat_to_netcdf()
c     -----------------------------------------------------  
c     Assemble and write track statistics / final state
c     Variables:
c         is_out_of_domain  =  0(inside) OR 1(outofdomain)
c         is_settled        = -1(dead)  OR 0(free) OR int>0(settled in box number)
c         degree_days       = degree-day integral for this particle
c         min_temp          = min temperature encountered for this particle
c         max_temp          = max temperature encountered for this particle
c         min_salt          = min salinity encountered for this particle
c         max_salt          = max salinity encountered for this particle
c     ----------------------------------------------------- 
      integer             :: i, isout,sett
c     -----------------------------------------------------
      do i=1,par_ens%last
c        ------------         
         if (par_ens%space_stack(i)%outofdomain) then
            isout = 1
         else
            isout = 0
         endif
         call nfcheck( nf90_put_var(ncid, is_out_of_domain_id,
     +        isout, start=(/i/)) )
c        ------------               
         if     (par_ens%state_stack(i)%type == particle_free) then
            sett = 0
         elseif (par_ens%state_stack(i)%type == particle_settled) then
            sett = par_ens%state_stack(i)%settleBox
         elseif (par_ens%state_stack(i)%type == particle_dead) then
            sett = -1
         else
            write(*,*) "save_tracer_stat_to_netcdf: unmapped state"
            stop 21
         endif   
         call nfcheck( nf90_put_var(ncid, is_settled_id,
     +        sett, start=(/i/)) )
c        ------------
         call nfcheck( nf90_put_var(ncid, age_at_settling_id,
     +        par_ens%state_stack(i)%bio%age, start=(/i/)))
c     ------------
         call nfcheck( nf90_put_var(ncid, degree_days_id,
     +        par_ens%state_stack(i)%bio%degree_days, start=(/i/)))
c        ------------
         call nfcheck( nf90_put_var(ncid, min_temp_id,
     +        par_ens%state_stack(i)%bio%min_temp, start=(/i/)))
c        ------------
         call nfcheck( nf90_put_var(ncid, max_temp_id,
     +        par_ens%state_stack(i)%bio%max_temp, start=(/i/)))       
c        ------------
         call nfcheck( nf90_put_var(ncid, min_salt_id,
     +        par_ens%state_stack(i)%bio%min_salt, start=(/i/)))
c        ------------
         call nfcheck( nf90_put_var(ncid, max_salt_id,
     +        par_ens%state_stack(i)%bio%max_salt, start=(/i/)))                
      enddo
      end subroutine save_tracer_stat_to_netcdf
      
      end program
