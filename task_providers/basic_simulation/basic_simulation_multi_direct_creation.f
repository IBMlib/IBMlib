ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Standard particle tracking simulation 
c     ---------------------------------------------------
c
c     setup and run an ensemble where particles are created by direct operator make_particles_direct
c     based on basic_simulation_direct_creation.f
c    
c     $Rev$
c     $LastChangedDate$
c     $LastChangedBy$ 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program tracker
      use input_parser
      use time_tools
      use run_context
      use physical_fields
      use particles
      use write_netcdf_lonlat_data
      
      implicit none
c     ------------ declarations ------------      
      integer      :: idum4(4)
      type(clock),target  :: start_time, end_time
      type(clock),pointer :: current_time
      type(particle_ensemble)     :: par_ens
      integer :: year, month, day, julday,last
      integer :: i,ix,iy,iz,dum(79),ilast,iunit,ihit
      integer :: idum(4), timeout, istep, srclabel
      integer :: isrc,nsrc,i2(2)
      integer :: start(500),nwords ! 500 arbitrary
      real    :: time_step, now,r4(4), pos(3)
      real    :: amp,xcen,ycen,vdiff,hdiff
      character*256 :: filename
      character*999 :: strbuf
      character(len=64) :: vattr(2,2) 
      type(netcdf_grid_2Ddata_handler) :: nchout
      
      integer :: nxs, nys, noutside
      real    :: lon1,lat1, dlon, dlat
      logical :: inside
      real, allocatable :: litter(:,:), area(:,:),xyz(:,:)
      integer, allocatable :: npar(:)
      character*999, allocatable :: initialization_data(:)
      
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
      call update_physical_fields()  

c     ---- create all particles once before main time loop ---
c
c     raw particle creation is at a fairly low level because many 
c     variants are anticipated; for simplicity just read one creation site
c
      nsrc  = count_tags(simulation_file, "create_particles")
      if (nsrc < 1) then
         write(*,*) "no particles - no run"
         stop 129
      else
         write(*,*) "found", nsrc, "particle sources"
      endif
      allocate( xyz(3,nsrc) )  
      allocate( npar(nsrc)  )
      allocate( initialization_data(nsrc)  )
c      
c     ------- read and parse all create_particles entries -------
c
      ihit    = 1
      do isrc = 1, nsrc
         call read_control_data(simulation_file, "create_particles",
     +                          strbuf, ihit)
         call tokenize(strbuf, start, nwords)
         if      (nwords < 4)  then
            
            stop "usage: create_particles = lon,lat,depth,nparticles"
         else if (nwords >= 4) then 
            read(strbuf,*) xyz(:,isrc), npar(isrc)
            if (nwords == 4) then
               initialization_data(isrc) = "" 
            else
               initialization_data(isrc) = strbuf(start(5):) ! rest of line contains bio init data
            endif
         endif
         write(*,743) isrc, xyz(1:2, isrc), npar(isrc)
         ihit = ihit + 1   ! take next box in control file
      enddo                     ! isrc = 1, nsrc
 743  format("source",i," @",1x,2f9.4," emits",i,1x,"particles")
c
c     ------- create particles (source label according to appearence in input) -------
c
      call setup_ensemble(par_ens, sum(npar)) 
      do isrc = 1, nsrc
         
         call generate_particles(par_ens, xyz(:,isrc), time_step,
     +        npar(isrc), isrc, initialization_data(isrc))
         
      enddo  ! isrc = 1, nsrc
c
c     ------- read output options -------
c      
      call read_control_data(simulation_file, "output_file", filename)
      call read_control_data(simulation_file, "output_grid_dims",i2)    ! nxs,nys
      nxs = i2(1)
      nys = i2(2)
      call read_control_data(simulation_file, "output_grid_bbox",r4)    ! W/S/E/N 
      dlon = (r4(3)-r4(1))/i2(1)
      dlat = (r4(4)-r4(2))/i2(2)
      lon1 = r4(1) + 0.5*dlon    ! cell center of first grid point
      lat1 = r4(2) + 0.5*dlat    ! cell center of first grid point 
              
      
      current_time => get_master_clock()
c     =====================  main time loop =====================
      istep   = 0
      
      do while (compare_clocks(current_time, end_time)
     +                         *nint(time_step) <= 0)             ! opposite for forward/backward simulation

         write(*,372) istep
         call update_physical_fields()
c        -------- propagate tracers  --------        
         call update_particles(par_ens, time_step)
         
c        -------- write tracer state -------- 
c        -------- loop control       --------
         call add_seconds_to_clock(current_time, nint(time_step))
         istep = istep + 1

      enddo
c      
c     -------- dump final simulation output as netcdf distribution --------
c
      vattr(1,1) = "long_name"
      vattr(2,1) = "probability density"
      vattr(1,2) = "unit"
      vattr(2,2) = "km^(-2)"
      call init_netcdf_dataset_write(nchout, filename,  
     +     nxs, nys, lon1, lat1, dlon, dlat,
     +     "xy", vname="litter_density",vattr=vattr)
     
      allocate ( litter(nxs, nys) )
      allocate ( area(nxs, nys)   )
      litter = 0
      noutside = 0
      call get_last_particle_number(par_ens, last)
      do i=1,last
         call get_particle_position(get_particle(par_ens,i), pos)
         call get_subgrid_indices(nchout, pos, ix, iy, inside)
         if (inside) then
            litter(ix, iy) = litter(ix, iy) + 1.0
         else
            noutside = noutside + 1
         endif
      enddo
      write(*,*) "fraction of particles outside analysis grid =",
     +     1.0*noutside/last

c     ------ write final data ------      
      call get_area_matrix(nchout, area)
      litter = litter/sum(litter)/area
      call write_netcdf_data(nchout, litter)
      call close_netcdf_dataset(nchout)
 372  format(25("*"), " main time loop step ",i5.5, " ", 25("*"))
 
c     ----------------- close down ---------------------------
     
      call close_physical_fields()
      call close_particles()
      deallocate ( litter )
      deallocate ( area )
      
      write(*,*) "normal termination of simulation"

      end program
