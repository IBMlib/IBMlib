ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ===============================================================
c            Linear Flow Fields Oceanography Provider
c     --------------------------------------------------------------- 
c     $Revision$
c     $LastChangedDate$
c     $LastChangedBy$ 
c     --------------------------------------------------------------- 
c     Provides a basic set of linear flow fields which can be used
c     for the developement and testing of code, independent of variable
c     oceanographic conditions. Flow fields have a set of default configs
c     that can be overrulled by options supplied via the controlfile, via
c     which it is possible to setup any combination of flat field, 
c     shear or helical flow.
c
c     Fields are defined by specifying two velocity vectors (a speed and
c     a direction), one at the top of the water column, and one at 
c     the bottom. A depth of the water column should also be specified
c     If the bottom vector is not specified, it defaults to the value
c     at the surface
c
c     Oceanography provider corresponds to a wet sphere - there is no
c     land, although the water column is finite.
c     =============================================================== 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module physical_fields
      use run_context, only: simulation_file
      use input_parser
      use time_tools           ! import clock type
      use constants

      implicit none
      private     

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
      public :: horizontal_range_check
      public :: coast_line_intersection

      interface  read_control_robust
        module procedure read_control_robust_vec
        module procedure read_control_robust_scalar
      end interface

c     -------------------- module data --------------------  
c
      type(clock), target  :: master_clock ! simulation master clock

c     -------------------- flow field variables -----------  
c
      real  :: surface_vel(2)  !speed (m/s), angle (degrees clockwise from North)
      real  :: bottom_vel(2)   !speed (m/s), angle (degrees clockwise from North)
      real  :: depth        ! m
      real  :: temp         ! deg C - uniform throughout 
      real  :: salin        ! psu   - uniform throughout
      real  :: hdiffus      ! horizontal diffusivity m^2/s
      real  :: vdiffus      ! vertical diffusivity   m^2/s
      
c     --------
      contains
c     --------

      subroutine init_physical_fields(time)
c     ---------------------------------------------------
c     module initialization
c     ---------------------------------------------------
      type(clock), intent(in),optional :: time
c     ---------------------------------------------------

c.....Display version numbers for reference
      write(*,*) "Linear Flow Fields Oceanography Provider: " //
     &   "$Rev$"
     
c.....Import configuration from control file, or default if not specified
      call read_control_robust(simulation_file,"surface_vel",
     &   surface_vel,(/1.0,90.0/))
      call read_control_robust(simulation_file,"bottom_vel",
     &   bottom_vel,surface_vel)   !Bottom defaults to the surface
      call read_control_robust(simulation_file,"depth",depth,1000.0)
      call read_control_robust(simulation_file,"temp",temp,15.0)
      call read_control_robust(simulation_file,"salin",salin,35.0)
      call read_control_robust(simulation_file,"hdiffus",hdiffus,10.0)
      call read_control_robust(simulation_file,"vdiffus",vdiffus,0.01)

c.....Report results
      write(*,*) "Surface velocity vector (m/s,deg) : ",surface_vel
      write(*,*) "Bottom velocity vector (m/s, deg) : ", bottom_vel
      write(*,*) "Water column depth (m)            : ", depth
      write(*,*) "Water temperature (deg C)         : ", temp
      write(*,*) "Salinity (psu)                    : ", salin
      write(*,*) "Horizontal diffusivity (m²/s)     : ", hdiffus
      write(*,*) "Vertical diffusivity (m²/s)       : ", vdiffus      
      
c.....Set clock if present     
      if (present(time)) master_clock = time
      end subroutine 

      subroutine close_physical_fields()
      end subroutine 


      function   get_master_clock()
      type(clock), pointer :: get_master_clock
      get_master_clock => master_clock
      end function 


      subroutine set_master_clock(time)
      type(clock), intent(in),optional :: time
      master_clock = time
      end subroutine 
  

      subroutine update_physical_fields(time, dt)
      type(clock), intent(in),optional :: time
      integer, intent(in),optional     :: dt
      end subroutine 


      subroutine interpolate_turbulence(xyz, r3, status)
c     ---------------------------------------------------
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r3(:)
      integer, intent(out) :: status
c     ---------------------------------------------------
      status =0
      r3=(/hdiffus,hdiffus,vdiffus/)
      end subroutine 


      subroutine interpolate_turbulence_deriv(xyz, r3, status)
c     ---------------------------------------------------
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r3(:)
      integer, intent(out) :: status
c     ---------------------------------------------------
      status =0
      r3=0.0
      end subroutine 

      subroutine interpolate_currents(xyz, r3, status)
c     ---------------------------------------------------
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r3(:)
      integer, intent(out) :: status
c     ------locals-------
      real ::   z, angle, speed,vel(2)
c     ---------------------------------------------------
      status =0
      z=xyz(3)
c     Interpolate velocity vector at depth
      vel=z/depth*(bottom_vel-surface_vel) + surface_vel
      speed=vel(1)
      angle =vel(2)
      r3(1)=speed*sin(deg2rad*angle)
      r3(2)=speed*cos(deg2rad*angle)
      r3(3)=0
      end subroutine 

      subroutine interpolate_temp (xyz, r, status) 
c     ---------------------------------------------------
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
c     ---------------------------------------------------
      status =0
      r=temp
      end subroutine 

      subroutine interpolate_salty(xyz, r, status)
c     ---------------------------------------------------
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
c     ---------------------------------------------------
      status =0
      r=salin
      end subroutine 

      subroutine interpolate_wind(xy, r2, status)
c     ---------------------------------------------------
      real, intent(in)     :: xy(:)
      real, intent(out)    :: r2(:)
      integer, intent(out) :: status
c     ---------------------------------------------------
      status =0
      r2=0.0
      end subroutine 

      subroutine interpolate_wdepth(xy, r, status) 
c     ---------------------------------------------------
      real, intent(in)     :: xy(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
c     ---------------------------------------------------
      status =0
      r=depth
      end subroutine 


      LOGICAL function is_wet(xyz)
      real, intent(in)     :: xyz(:)
      if(xyz(3) < depth) then
        is_wet = .TRUE.
      else
        is_wet = .FALSE.
      endif
      end function

      LOGICAL function is_land(xy)
      real, intent(in)     :: xy(:)
      is_land = .FALSE.
      end function

      LOGICAL function horizontal_range_check(xy)
      real, intent(in)     :: xy(:)
      horizontal_range_check = .TRUE.
      end function

      subroutine coast_line_intersection(xyz1, xyz2, anycross, xyzref, 
     +                                   xyzhit) 
c     ---------------------------------------------------
c     xyzref and xyzhit do not need to be defined, as there
c     is no coastline here - anycross will always be false
c     ---------------------------------------------------
      real, intent(in)     :: xyz1(:),xyz2(:)
      logical, intent(out) :: anycross
      real, intent(out)    :: xyzref(:), xyzhit(:)
c     ---------------------------------------------------
      anycross = .FALSE.
      xyzref =0.0
      xyzhit = 0.0
      end subroutine 


      subroutine read_control_robust_vec(ctrlfile,tag,variable,def)
c     ---------------------------------------------------
c     wrapper to the standard read_control_data
c     that is robust to the absence of values in the
c     controlfile - in that case, the value returned is
c     the default argument
c     ---------------------------------------------------
      character(*)        :: tag
      type(control_file)   :: ctrlfile
      real                 :: variable(:),def(:)
c     -------locals--------
      integer  ntags
c     ---------------------------------------------------
      ntags = count_tags(ctrlfile,tag)
      if(ntags>0) then
        call read_control_data(ctrlfile,tag,variable)
      else
        variable=def
      endif
      end subroutine

      subroutine read_control_robust_scalar(ctrlfile,tag,variable,def)
c     ---------------------------------------------------
c     wrapper to the standard read_control_data
c     that is robust to the absence of values in the
c     controlfile - in that case, the value returned is
c     the default argument
c     ---------------------------------------------------
      character(*)        :: tag
      type(control_file)   :: ctrlfile
      real                 :: variable,def
c     -------locals--------
      integer  ntags
c     ---------------------------------------------------
      ntags = count_tags(ctrlfile,tag)
      if(ntags>0) then
        call read_control_data(ctrlfile,tag,variable)
      else
        variable=def
      endif
      end subroutine
      
      end module
