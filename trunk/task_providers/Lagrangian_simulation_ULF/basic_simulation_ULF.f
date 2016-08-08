cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     -------------------------------------------------------------------------
c     Standard particle tracking simulation 
c     writing output in Unified Lagrangian Format (ULF)
c     as stored in netcdf :
c
c       section 0) premises (input, version, configurations, ...) - not implemented yet
c       section 1) optional topography (for visualization)
c       section 2) optional particle tracks (number of frames bind to unlimited dimension)
c       section 3) optional connectivity matrix                   - not implemented yet
c
c     This version implement sections 1+2
c     -------------------------------------------------------------------------
c     https://www.unidata.ucar.edu/software/netcdf/netcdf-4/newdocs/netcdf-f90/
c     -------------------------------------------------------------------------
c     $Rev: 175 $
c     $LastChangedDate: 2011-01-06 17:23:41 +0100 (Thu, 06 Jan 2011) $
c     $LastChangedBy: asch $ 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program tracker
      use netcdf
      use input_parser
      use time_tools
      use run_context
      use physical_fields
      use particles
     
      implicit none
c     ------------ declarations ------------      
      type(clock),target  :: start_time, end_time
      type(clock),pointer :: current_time
      type(particle_ensemble)     :: par_ens
      type(emission_box), pointer :: emitboxes(:)
      integer :: year, month, day, julday
      integer :: i,ix,iy,iz,dum(79),ilast,iunit, numsec, npart
      integer :: idum0(4), idum1(4), timeout, istep, ntracers, ntim
      integer :: ncid, nskip, topolon_id, topolat_id, npart_id
      integer :: lon_id, lat_id, depth_id, time_id, i4_id, nframes_id
      integer :: nxtop, nytop, nxtop_id, nytop_id, topography_id
      real    :: time_step, t, rdum6(6), dt_save
      real, allocatable :: topography(:,:), topolon(:), topolat(:)
      logical :: save_tracks, save_xy, save_xyz, topo_scan
      character*256 :: filename, topofilename
      

c     ------ for windows debugging bind stdout/stderr to specific files (non portable) ------

      open(0, file="stderr.txt")
      open(6, file="stdout.txt")

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
      call update_physical_fields()  ! need topology to verify emission boxes
      call setup_ensemble_from_file(par_ens, "emitbox", emitboxes)
      npart = get_ensemble_size(par_ens) ! at this point we know max number of particles
c
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
      else  ! section 1 empty
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
         call nfcheck( nf90_put_att(ncid, NF90_GLOBAL, 
     +                 "topography_file", topofilename) )
      elseif (count_tags(simulation_file, "save_xy_frequency") > 0) then
         call read_control_data(simulation_file, "save_xy_frequency", 
     +                     dt_save)
         save_xy  = .true.
         save_xyz = .false.
      else
         save_xy  = .false.
         save_xyz = .false.
         dt_save  = 0.0  ! render defined, used below
      endif
      save_tracks = save_xy .or. save_xyz
      nskip = max(1, nint(dt_save/time_step))  ! skip frames between saves

c     NF90_UNLIMITED binds to particle frames

      if (save_tracks) then
         call nfcheck( nf90_def_dim(ncid, "nframes", NF90_UNLIMITED, 
     +                              nframes_id))     
         call nfcheck( nf90_def_dim(ncid, "nparticles",npart, npart_id))
         call nfcheck( nf90_def_dim(ncid, "i4",        4,     i4_id))
c        -------- define lon --------
         call nfcheck( nf90_def_var(ncid, "lon", NF90_FLOAT, 
     +                 (/npart_id, nframes_id/), lon_id))
         call nfcheck( nf90_put_att(ncid, lon_id, "meaning", 
     +                 "particles zonal position")) 
         call nfcheck( nf90_put_att(ncid, lon_id, "unit", 
     +                 "degrees East")) 
c        -------- define lat --------
         call nfcheck( nf90_def_var(ncid, "lat", NF90_FLOAT, 
     +                 (/npart_id, nframes_id/), lat_id))
         call nfcheck( nf90_put_att(ncid, lat_id, "meaning", 
     +                 "particles meridional position")) 
         call nfcheck( nf90_put_att(ncid, lat_id, "unit", 
     +                 "degrees North")) 
c        -------- define time --------
         call nfcheck( nf90_def_var(ncid, "time", NF90_INT, 
     +                 (/i4_id, nframes_id/), time_id))
         call nfcheck( nf90_put_att(ncid, time_id, "meaning", 
     +                 "time corresponding to particles frames")) 
         call nfcheck( nf90_put_att(ncid, time_id, "unit", 
     +                            "year, month, day, second_in_day"))
      endif
      if (save_xyz) then
c        -------- define depth --------
         call nfcheck( nf90_def_var(ncid, "depth", NF90_FLOAT, 
     +        (/npart_id, nframes_id/), depth_id))   
         call nfcheck( nf90_put_att(ncid, depth_id, "meaning", 
     +                 "particles vertical position")) 
         call nfcheck( nf90_put_att(ncid, depth_id, "unit", 
     +                              "meters below sea surface"))
      endif

c
c     --- leave define mode - after here we can write data variables
c
      call nfcheck( nf90_enddef(ncid) )  
      
      if (topo_scan) then
         call scan_topography() ! topolon,topolat -> topography
         call nfcheck( nf90_put_var(ncid, topolon_id, topolon) )
         call nfcheck( nf90_put_var(ncid, topolat_id, topolat) )
         call nfcheck( nf90_put_var(ncid, topography_id, topography) )
      endif
      
      

      current_time => get_master_clock()
c     =====================  main time loop =====================
      t = 0.0
      do istep = 1, ntim
         write(*,372) istep
         if ( save_tracks .and. (mod(istep-1, nskip) == 0)) then
            call write_particle_frame()
         endif
         call update_physical_fields()      
         call generate_particles(par_ens, emitboxes, time_step)
         call update_particles(par_ens, time_step)
         t  = t + time_step
         call add_seconds_to_clock(current_time, nint(time_step))
      enddo

    
 372  format(25("*"), " main time loop step ",i5.5, " ", 25("*"))

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

      end program
