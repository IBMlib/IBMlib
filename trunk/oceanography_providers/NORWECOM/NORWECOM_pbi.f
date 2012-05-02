ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     -----------------------------------------------------------
c     NORWEgian ECOlogical Model (NORWECOM) Oceanography Provider
c     -----------------------------------------------------------
c     $Rev$
c     $LastChangedDate$
c     $LastChangedBy$ 
c
c     This Oceanography Provider uses the generic mesh_grid pbi element
c     to provide a general interface to the NORWECOM data set.
c
c    TODO:--------------------------------
c       * Check calculation of horizontal diffusivity
c       * Check coastline reflection
c    TODO later:--------------------------
c       * Check that velocity interpolation is appropriate
c         especially around the edge
c       * Test velocity vectors by comparison with ARIANE
c       * Improve interpolations around coastline
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module physical_fields
      use mesh_grid,
     +   ncc_interpolate_currents  => interpolate_currents,
     +   ncc_interpolate_turbulence=> interpolate_turbulence,
     +   ncc_interp_turb_deriv     => interpolate_turbulence_deriv
      use time_services
      use horizontal_grid_transformations
      use run_context, only: simulation_file
      use netcdf
      use array_tools
      use input_parser
      use constants
      implicit none
      private                     ! default visibility

c     -------------------- public interface --------------------  
      public :: init_physical_fields
      public :: close_physical_fields
      public :: get_master_clock
      public :: set_master_clock
      public :: update_physical_fields
      public :: interpolate_turbulence
      public :: interpolate_turbulence_deriv
      public :: interpolate_currents
      public :: interpolate_temp 
      public :: interpolate_salty
      public :: interpolate_wind
      public :: interpolate_wdepth
      public :: is_wet    
      public :: is_land
      public :: coast_line_intersection
      public :: horizontal_range_check
      public :: get_pbi_version

c     -------------------- helper arrays ------------------  
      logical,allocatable :: land(:,:)
      real,allocatable    :: sigma_edges(:)
      real,allocatable    :: sigma_lwidth(:)

c     -------------------- module data --------------------  
      character*999 :: hydroDBpath ! hydrographic data sets
      character*1, parameter ::  path_separator = "/"   ! UNIX/Linux separator
      real, parameter :: NaN = 0.0 / 0.0
      integer :: verbose = 0 !0=all off. 1=updates only. 2=all

c     Time pieces
c     We are forced to make an ugly hack here, as our time library can't 
c     handle times before 1970, but the NORWECOM baseline is set to 1 Jan 
c     1950. We therefore need to have a constant integer that accounts for 
c     the difference between these two, NORWECOM_t_offset. The value of this
c     integer is the time, in POSIX time, where NORWECOM time units are based
      integer, parameter     ::  NORWECOM_t_offset = 631148400 
      type(clock), private, target :: master_clock
      type(clock)  ::  cur_daily_fields, cur_hourly_fields 


c     ===================================================
                         contains
c     ===================================================
      include "NORWECOM_interpolators.f"

      character*200 function get_pbi_version()  
      get_pbi_version = "NORWECOM pbi. $Rev$ " 
      end function


      subroutine init_physical_fields(time)
c     ---------------------------------------------------
c     Initialises the physical fields module, including
c     allocating the data arrays
c     ---------------------------------------------------
      character*999        :: bath_fname
      integer              :: ncid, varid,ix,iy,iz,topo_fill
      integer,allocatable  :: intbath(:,:)
      type(clock), intent(in),optional :: time
c     ---------------------------------------------------
c.....Display version numbers for reference
      write(*,*) trim(get_pbi_version())

c.....Setup the grid. 
      nz=20
      call init_horiz_grid_transf()

c.....Initialise the mesh grid handler
      call init_mesh_grid(NaN)

c.....Set clock if present     
      if (present(time)) call set_master_clock(time)

c.....Read input path      
      call read_control_data(simulation_file,"hydroDBpath",hydroDBpath)
      write(*,*) "init_physical_fields: hydrographic database path =", 
     +           trim(hydroDBpath)
 

c.....Read data from the bathmetry file
      allocate(intbath(nx,ny)) !Temporary integer bathymetry 
      allocate(sigma_edges(nz+1))
      allocate(sigma_lwidth(nz))
      call read_control_data(simulation_file,"bath_file", bath_fname) 
      write(*,*) "init_physical_fields: bathymetry file = ",
     +    trim(adjustl(bath_fname))
      call check_ncdf(nf90_open(trim(adjustl(bath_fname)),
     +        nf90_nowrite,ncid))
      call check_ncdf(nf90_inq_varid(ncid,"Topo",varid))
      call check_ncdf(nf90_get_var(ncid,varid,intbath))
      call check_ncdf(nf90_get_att(ncid,varid,"_FillValue",topo_fill))
      call check_ncdf(nf90_inq_varid(ncid,"Zedges",varid))
      call check_ncdf(nf90_get_var(ncid,varid,sigma_edges))
      call check_ncdf(nf90_close(ncid))

c.....setup bathymetry query arrays                
c     The NORWECOM grid is terrain following, so every point is wet, except those on land.
      allocate(land(nx,ny))    !Array containing true/false for land
      bottom_layer =nz
      wetmask=1   !1 is wet, 0 is dry
      land = .FALSE.
      wdepth = real(intbath)
      where(intbath==topo_fill )!Bathymetry is fill value for dry land. 
          wetmask=0         !Mark as dry 
          land = .TRUE.
          bottom_layer=0    !bottom layer is the 0th layer
          wdepth = NaN      !Set depth to 0
      end where
      deallocate(intbath)

c.....setup vertical grid auxillary arrays
      acc_width(:,:,1)=0.
      do iz=1,nz
         sigma_lwidth(iz) = sigma_edges(iz+1)-sigma_edges(iz)
         ccdepth(:,:,iz)=wdepth*(sigma_edges(iz+1)+sigma_edges(iz))*0.5
         acc_width(:,:,iz+1)=wdepth*sigma_edges(iz+1)
      enddo 
      end subroutine init_physical_fields
      
      
      subroutine close_physical_fields()
c-------------------------------------------------------
c     Cleanup module
c-------------------------------------------------------
      deallocate(land)
      deallocate(sigma_edges)
      deallocate(sigma_lwidth)
      end subroutine close_physical_fields     

      

      subroutine update_physical_fields(force)
c     -----------------------------------------------------
c     Update physical_fields corresponding to match master clock
c     The decision to update  is based upon the distance between the 
c     current time and the current loaded fields - if this is greater 
c     than the time resolution  of the corresponding  fields, an 
c     update is required. An update can be forced by setting the
c     force parameter to TRUE
c     -----------------------------------------------------
      logical,  optional :: force
      integer       :: cyear, cmonth, cday, cisec
      integer       :: offset
      logical       :: updated
      type(clock), pointer :: ctime
      type(timer)   :: stopwatch
c     -----------------------------------------------------
      !Start stopwatch to time update
      call tic(stopwatch) 

      !Set pointer to clock
      ctime => get_master_clock()

c     First if we are forcing the update, lets jump straight there
      if(present(force))  then
      if(force) then
        call logmsg(2,"Forcing update...")
        call update_daily_physical_fields(ctime)
        call update_hourly_physical_fields(ctime)
        call update_derived_variables()
        if(verbose>=1) then
          call toc(stopwatch,"Physical_fields update complete :") 
        endif
        return
      endif
      endif 

c.....Decide whether an update of the daily fields is needed or not
      updated=.FALSE.
      call get_period_length_sec(ctime,cur_daily_fields,offset)
      if(abs(offset) > 43200)  then
        call logmsg(2,"Updating daily fields....")
        call update_daily_physical_fields(ctime)
        updated = .TRUE.
      endif

c.....Now decide whether an update of the hourly fields is needed or not,
c     using similar logic
      call get_period_length_sec(ctime,cur_hourly_fields,offset)
      if(abs(offset) > 1800 ) then
        call logmsg(2,"Updating hourly fields....")
        call update_hourly_physical_fields(ctime)
        updated = .TRUE.
      endif

c.....If fields were updated then update dependent variables      
      if (updated) then
        call update_derived_variables()
        !How long did it all take?
        if(verbose>=1) then
          call toc(stopwatch,"Physical_fields update complete :") 
        endif
      endif
      
      end subroutine update_physical_fields


      subroutine update_daily_physical_fields(ctime)
c     -----------------------------------------------------
c     Update physical_fields corresponding to clock time ctime
c     We assume that the NORWECOM daily fields represent the
c     entirity of the given day. Thus, the first entry, which
c     is timestamped with 2006/01/01 12:00 covers everything
c     from 00:00 to 23:59
c     -----------------------------------------------------
      type(clock), intent(in)   :: ctime
      type(clock)   :: ctmp
      integer       :: cyear, cmonth, cday
      integer       :: jday,idx,leng
      real          :: actual_idx(1) !Index that we actually got
c     -----------------------------------------------------
c.....Figure out what the next data frame should be
      call get_date_from_clock(ctime,cyear,cmonth,cday)
      call get_julian_day(ctime,jday)  !The julian day is the index used to store the daily fields
      idx=jday    ! Index used to look up data

c.....Now load the new daily fields
      call read_NORWECOM_fields("STbio_daily","Temp",cyear,idx,temp)
      call land_fill(temp,NaN)  !Set land to NaN
      call read_NORWECOM_fields("STbio_daily","Salt",cyear,idx,salinity)
      call land_fill(salinity,NaN) !Set land to NaN

c.....Get the NORWECOM time that the array index corresponds to      
      call read_NORWECOM_fields("STbio_daily","T",cyear,idx,
     +         dat1d=actual_idx)

c.....And finish off the update by setting the time properly etc and
c     checking it agrees with result from NetCDF files
      call set_clock(ctmp,cyear,cmonth,cday,43200)
      leng = get_POSIXtime(ctmp)
      if (nint(actual_idx(1))*3600.ne.leng+NORWECOM_t_offset) then
          write(*,'(25("-"),a,25("-"))') "STOP"
          write(*,*)"Update_daily_physical_ fields: data frame" //
     +       " loaded does not match that desired."
          write(*,*) "Actual time : ",get_datetime(ctime)
          write(*,*) "Desired time frame : ",get_datetime(ctmp)
          write(*,*) "Index of loaded time frame : ",actual_idx  
          stop
      else
          call set_clock_from_clock(cur_daily_fields,ctmp)
      endif
c     Confirm update
      call logmsg(2, "New daily fields : "
     &                  //get_datetime(cur_daily_fields))
 
      end subroutine update_daily_physical_fields


      subroutine update_hourly_physical_fields(ctime)
c     -----------------------------------------------------
c     Update physical_fields corresponding to clock time ctime
c     We assume that the NORWECOM fields are representative of 
c     the hour period in which they lie in the middle of  
c     eg the first field in the NCDF file,
c     which is time stamped 2006/01/01 01:00 represents the flow
c     fields between 00:30 and 01:30.
c     -----------------------------------------------------
      type(clock), intent(in)   :: ctime
      type(clock)   :: ctmp, cstart_of_year
      integer       :: cyear, cmonth, cday,csec
      integer       :: jday, idx, secon
      integer       :: i,j,leng, elapsed
      real          :: actual_idx(1)
c     -----------------------------------------------------
c.....Figure out what the next data frame should be
      call set_clock_from_clock(ctmp,ctime)
      call add_seconds_to_clock(ctmp,-1800) !Taking the clock back some 30min avoids problems at the very start of the year
      call get_date_from_clock(ctmp,cyear,cmonth,cday)
      call set_clock(cstart_of_year,cyear,1,1,0)  
      call get_period_length_sec(ctmp,cstart_of_year,elapsed) !Time elapsed since the start of the year (in secs)
      idx = ceiling(real(elapsed)/3600.)
      if (idx.eq.0) then 
        idx=1
      end if

c.....Now load the new hourly fields
      call read_NORWECOM_fields("velo_hourly","U",cyear,idx,u)
      call land_fill(u,0.0)
      call read_NORWECOM_fields("velo_hourly","V",cyear,idx,v)
      call land_fill(v,0.0)
      call read_NORWECOM_fields("velo_hourly","VertD",cyear,
     +        idx,vdiffus(:,:,1:nz))

c.....Polish up the vertical diffusivity - diffusivity at the boundaries
c     conditions should be zero. Otherwise, diffusivity should not
c     fall below the minimum background value of 1e-5.
      where (vdiffus(:,:,2:nz) < 1e-5)
          vdiffus(:,:,2:nz) = 1e-5
      end where
      vdiffus(:,:,1) =0
      vdiffus(:,:,nz+1) =0
      call land_fill(vdiffus(:,:,1:nz),0.0)

c.....Get the NORWECOM time that the array index corresponds to      
      call read_NORWECOM_fields("velo_hourly","T",cyear,idx,
     +         dat1d=actual_idx)

c.....And finish off the update by setting the time properly etc and
c     checking it agrees with result from NCDF files
      call set_clock_from_clock(ctmp,cstart_of_year)
      call add_seconds_to_clock(ctmp,idx*3600)
      leng= get_POSIXtime(ctmp)
      if (nint(actual_idx(1))*3600.ne.leng+NORWECOM_t_offset) then
          write(*,'(25("-"),a,25("-"))') "STOP"
          write(*,*)"Update_hourly_physical_ fields: data frame" //
     +       " loaded does not match that desired."
          write(*,*) "Actual time : ", get_datetime(ctime)
          write(*,*) "Desired time frame : ",get_datetime(ctmp)
          write(*,*) "Index of loaded time frame : ",actual_idx  
          stop
      else
          call set_clock_from_clock(cur_hourly_fields,ctmp)
      endif
c     Confirm update
      call logmsg(2,"New hourly fields : "
     +             //get_datetime(cur_hourly_fields))

      end subroutine update_hourly_physical_fields


      subroutine read_NORWECOM_fields(prefix,varName,year,idx,
     +                  dat3d,dat2d,dat1d)
c     -----------------------------------------------------
c     Reads data directly from the NORWECOM NetCDF database. The
c     subroutine accepts either a 3d, 2d or 1d output variable
c     which is then used to determine how the data is read
c     -----------------------------------------------------
      character(len=*),intent(in) ::prefix,varName
      integer,intent(in)          ::year,idx
      real,optional,intent(out)   ::dat3d(:,:,:)
      real,optional,intent(out)   ::dat2d(:,:)
      real,optional,intent(out)   ::dat1d(1)
      integer       :: ncid,varid
      character*999 :: ncfilename
      integer       :: status
      real          :: scf
c     -----------------------------------------------------
c     Setup filename
      write(ncfilename,801) trim(adjustl(hydroDBpath)),
     &        path_separator,prefix,year
 801  format(3a,i4,".nc")
      call logmsg(2,"Loading "//varName//" data from: "
     +      //trim(ncfilename))
c     Now open the file
      call check_ncdf(nf90_open(ncfilename,nf90_nowrite,ncid))
      call check_ncdf(nf90_inq_varid(ncid,varName,varid))
c     Get scale factor if it exists 
      scf= 1.0
      status = nf90_inquire_attribute(ncid,varid,"scale_factor")
      if (status.eq. nf90_noerr) then
        call check_ncdf(nf90_get_att(ncid,varid,"scale_factor",scf))
      end if
c     Read the data
      if (present(dat3d)) then 
         call check_ncdf(nf90_get_var(ncid,varid,dat3d,
     &            start=(/1,1,1,idx/),count=(/nx,ny,nz,1/) )) 
         dat3d=dat3d*scf
      elseif (present(dat2d)) then
         call check_ncdf(nf90_get_var(ncid,varid,dat2d,
     &            start=(/1,1,idx/),count=(/nx,ny,1/) )) 
         dat2d=dat2d*scf
      elseif (present(dat1d)) then
         call check_ncdf(nf90_get_var(ncid,varid,dat1d,
     &            start=(/idx/),count=(/1/) )) 
         dat1d=dat1d*scf
      else
         write(*,*) "Arguments to read_NORWECOM_fields misspecified"
         stop
      endif

c     Tidy up
      call check_ncdf(nf90_close(ncid))


      end subroutine read_NORWECOM_fields


      subroutine land_fill(dat,val)
c     -----------------------------------------------------
c     Sets the data stored in the supplied matrix to a 
c     value of "val" at all points that are on land
      real, intent(inout) :: dat(:,:,:)
      real, intent(in)    :: val
      integer ix, iy
c     -----------------------------------------------------
      do ix=1,nx
      do iy=1,ny
       if(wetmask(ix,iy)==0) dat(ix,iy,:) = val
      enddo
      enddo
      end subroutine land_fill


      subroutine check_ncdf(status)
c     -------------------------------------------------------
c     Checks that NetCDF operations worked
c     -------------------------------------------------------
      integer, intent ( in) :: status
    
      if(status /= nf90_noerr) then 
         print *, trim(nf90_strerror(status))
         call abort("Check_ncdf","Stopped")
      end if
      end subroutine check_ncdf


      subroutine update_derived_variables()
c     -----------------------------------------------------
c     Initiates a recalculation of derived variables the
c     flow fields have been updated
c     -----------------------------------------------------
      call logmsg(2,"Updating derived variables....")
      call update_w()
      call update_horizontal_turbulence()
      end subroutine update_derived_variables


      subroutine update_w()
c     -----------------------------------------------------
c     Update w(:,:,:) from u(:,:,:), v(:,:,:) using 
c     water incompressibility, starting at each bottom face
c     with w=0 and integrating upward
c     Notice that w(nx,:,:) and w(:,ny,:) can not be defined
c     by this procedure, due to open faces (these planes are 
c     assigned w=0)
c
c     w is midt point at cell upper face and positive upward 
c     so that sea bed boundary condition w=0 is implicit   
c     In principle, w(:,:,1) should equal (d/dt) dslm(:,:)
c     (assertion not checked)
c     -----------------------------------------------------
      integer :: ix, iy, iz
      real    :: dw, wlow, aoh, totwth
c     -----------------------------------------------------      
      w = 0.              ! pad value for undefined entries
      do ix=1,nx-1        ! w(nx,:,:)  undefined 
         do iy=1,ny-1     ! w(:,ny,:) undefined
            wlow = 0.     ! sea bed boundary condition
            do iz = bottom_layer(ix,iy), 1, -1  ! bottom_layer>= 1
               totwth = sigma_lwidth(iz)*wdepth(ix,iy)
               aoh = cell_area(ix,iy)/totwth
               dw  = yfacelen(ix,iy)  * u(ix,iy,iz) - 
     1               yfacelen(ix+1,iy)* u(ix+1,iy,iz) + 
     2               xfacelen(ix,iy)  * v(ix,iy,iz) - 
     3               xfacelen(ix,iy+1)* v(ix,iy+1,iz)
               w(ix,iy,iz) = wlow + dw/aoh
               wlow        = w(ix,iy,iz) ! is now next lower face
            enddo
         enddo
      enddo
  
      end subroutine update_w


      subroutine update_horizontal_turbulence()
c     -----------------------------------------------------
c     Smagorinsky parameterization of horizontal eddy diffusivity
c     Parameterisation employed here is based on that employed
c     natively in the NORWECOM code
c  
c     hdiffus is cell centered, i.e. hdiffus(ix,iy,iz) corresponds to (x,y,z) = (ix,iy,iz)
c     Evaluate shear by finite differencing (with one grid spacing from center position)
c     -----------------------------------------------------
c     Smagorinsky parameters for horizontal diffusion      
      real,parameter          :: smagorinsky_constant =0.03  !from M Skogen
c     local variables
      integer                 :: I,J,K, KB
c     ---------------------------------------------------------------------------
      call logmsg(2,"update_horizontal_turbulence: entered.")
c     Setup KB counter to match NORWECOM code
      KB = NZ +1  !Number of sigma  layers +1 
c     Populate the hdiffus array
      DO K = 1,KB-1
        DO J = 2,NY-1
          DO I = 2,NX-1
            hdiffus(I,J,K)=smagorinsky_constant*DX*DY
     1       *SQRT( ((U(I+1,J,K)-U(I,J,K))/DX)**2
     2       +((V(I,J+1,K)-V(I,J,K))/DY)**2
     3       +.5*(.25*(U(I,J+1,K)+U(I+1,J+1,K)-U(I,J-1,K)-U(I+1,J-1,K))
     4        /DY   
     5       + .25*(V(I+1,J,K)+V(I+1,J+1,K)-V(I-1,J,K)-V(I-1,J+1,K))
     6       /DX)**2)
          ENDDO   !I looping in horizontal x
        ENDDO     !J looping in horizotnal y
      ENDDO       !K looping in vertical direction z

      !Deviation from NORWECOM source code
      !Unit conversions. DX and DY in our setup are in km, but we
      !want hdiffus in m^/2. Unit analysis of the above expression shows it gives
      !hdiffus in units of km.m/s. Therefore multiple by 1000 m/km to get m^2/s
      hdiffus = hdiffus* 1000
      !/End deviation
C 
C     Extend the diffusivities to the boundary grid cells.
C
      DO  K = 1,KB-1
        DO  J = 2,NY-1
          hdiffus(1,J,K) = hdiffus(2,J,K)
          hdiffus(NX,J,K) = hdiffus(NX-1,J,K)
        ENDDO     !J looping in horizotnal y
      ENDDO       !K looping in vertical direction z
      DO  K = 1,KB-1
        DO  I = 2,NX-1
          hdiffus(I,1,K) = hdiffus(I,2,K)
          hdiffus(I,NY,K) = hdiffus(I,NY-1,K)
        ENDDO     !I looping in horizontal x
      ENDDO       !K looping in vertical direction z
      DO  K = 1,KB-1
        hdiffus(1,1,K) = hdiffus(2,2,K)
        hdiffus(1,NY,K) = hdiffus(2,NY-1,K)
        hdiffus(NX,NY,K) = hdiffus(NX-1,NY-1,K)
        hdiffus(NX,1,K) = hdiffus(NX-1,2,K)
      ENDDO       !K looping in vertical direction z
C
c     Original code has a 2D implementation here. We don't have anything
c     analogous in IBMlib, so ignore. However, have maintained code here
c     for consistency and reference
!!      DO 100 J=1,NY
!!        DO 100 I=1,NX
!!          AAM2D(I,J)=0.
!!100   CONTINUE
!!      DO 110 K = 1,KB-1
!!        DO 110 J = 1,NY
!!          DO 110 I = 1,NX
!!            AAM2D(I,J)=AAM2D(I,J)+AAM(I,J,K)*DZ(K)
!!110   CONTINUE
!!C
      end subroutine update_horizontal_turbulence


      subroutine logmsg(level,msg)
      integer, intent(in) :: level
      character(*), intent(in) :: msg
      if(verbose>=level) then
         write(*,*) trim(msg)
      endif
      end subroutine 



      end module
      

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$      PROGRAM test
c$$$C     ifort -e90 physical_fields.f run_context.o input_parser.o time_tools.o dmidata_interface.o libtime/libtime77.a gribex/libgribex.a
c$$$      use physical_fields
c$$$      integer :: i,j,k, ic,ix,iy,iz
c$$$      integer, parameter :: ns = 7
c$$$      real :: pos(3,ns**3), uvw(3,ns**3),x,y,z
c$$$c      real    :: spts(ns) = (/-0.7,  0.49, 0.51, 1.49,
c$$$c     +                         1.51, 2.49, 2.51, 8.0/)
c$$$  real    :: spts(ns) = (/1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1/)
c$$$
c$$$c     nx=2; ny=2; nz=2
c$$$      nx=4; ny=4; nz=4
c$$$
c$$$      allocate( u(nx,ny,nz)         ) 
c$$$      allocate( v(nx,ny,nz)         ) 
c$$$      allocate( w(nx,ny,nz)         ) 
c$$$      allocate( bottom_layer(nx,ny) ) 
c$$$
c$$$c$$$      u=0.
c$$$c$$$      v=0.
c$$$c$$$      w=0.
c$$$c     u=3x+1  v=4y+2  w=5z+7
c$$$      do ix=1,nx
c$$$         do iy=1,ny
c$$$            do iz=1,nz
c$$$               u(ix,iy,iz) = 3*(ix+0.5) + 1
c$$$               v(ix,iy,iz) = 4*(iy-0.5) + 2
c$$$               w(ix,iy,iz) = 5*(iz-0.5) + 7
c$$$            enddo
c$$$         enddo
c$$$      enddo
c$$$      bottom_layer=nz
c$$$c      call random_number(pos)
c$$$c      pos(1,:) = pos(1,:)*1.5*nx - 1.0
c$$$c      pos(2,:) = pos(2,:)*1.5*ny - 1.0
c$$$c      pos(3,:) = pos(3,:)*1.5*nz - 1.0
c$$$      ic=1
c$$$      do i=1,ns
c$$$         do j=1,ns
c$$$            do k=1,ns
c$$$               pos(:,ic) = (/spts(i), spts(j), spts(k)/)
c$$$               ic = ic+1
c$$$            enddo
c$$$         enddo
c$$$      enddo
c$$$      
c$$$      call interpolate_currents(pos,uvw)
c$$$
c$$$      do ic=1,ns**3
c$$$         x = pos(1,ic)
c$$$         y = pos(2,ic)
c$$$         z = pos(3,ic)
c$$$         write(*,33) x,y,z, uvw(1,ic)-(3*x+1), 
c$$$     +                     uvw(2,ic)-(4*y+2), uvw(3,ic)-(5*z+7)
c$$$ 33      format(6f8.3)
c$$$      enddo
c$$$
c$$$      end PROGRAM 


