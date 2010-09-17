ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Empty public interface of physical_fields
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module physical_fields
      use time_tools           ! import clock type
      implicit none

      interface get_local_distance
         module procedure get_local_distance_scalar
         module procedure get_local_distance_vector      
      end interface

      type(clock), private, target :: master_clock

      contains
      
    
      subroutine init_physical_fields(time)
      type(clock), intent(in),optional :: time
      end subroutine 

      subroutine close_physical_fields()
      end subroutine 

      function   get_master_clock()
      type(clock), pointer :: get_master_clock
      get_master_clock => master_clock
      end function 

      subroutine set_master_clock(time)
      type(clock), intent(in),optional :: time
      end subroutine 
  

      subroutine update_physical_fields(time, dt)
      type(clock), intent(in),optional :: time
      integer, intent(in),optional     :: dt
      end subroutine 

      subroutine interpolate_turbulence(xyz, r3, status)
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r3(:)
      integer, intent(out) :: status
      end subroutine 

      subroutine interpolate_turbulence_deriv(xyz, r3, status)
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r3(:)
      integer, intent(out) :: status
      end subroutine 

      subroutine interpolate_currents(xyz, r3, status)
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r3(:)
      integer, intent(out) :: status
      end subroutine 

      subroutine interpolate_temp (xyz, r, status) 
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
      end subroutine 

      subroutine interpolate_salty(xyz, r, status)
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
      end subroutine 

      subroutine interpolate_wind(xy, r2, status)
      real, intent(in)     :: xy(:)
      real, intent(out)    :: r2(:)
      integer, intent(out) :: status
      end subroutine 

      subroutine interpolate_wdepth(xy, r, status) 
      real, intent(in)     :: xy(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
      end subroutine 


      LOGICAL function is_wet(xyz)
      real, intent(in)     :: xyz(:)
      end function

      LOGICAL function is_land(xy)
      real, intent(in)     :: xy(:)
      end function

      LOGICAL function horizontal_range_check(xy)
      real, intent(in)     :: xy(:)
      end function

      subroutine get_horizontal_distance(xy1, xy2, r)
      real, intent(in)     :: xy1(:), xy2(:)
      real, intent(out)    :: r
      end subroutine 


      subroutine get_local_distance_scalar(xyz1, xyz2, r)
      real, intent(in)     :: xyz1(:),xyz2(:)
      real, intent(out)    :: r
      end subroutine 
      subroutine get_local_distance_vector(xyz1, xyz2, v)
      real, intent(in)     :: xyz1(:),xyz2(:)
      real, intent(out)    :: v(:)
      end subroutine 

      
      subroutine coast_line_intersection(xyz1, xyz2, anycross, xyzref, 
     +                                   xyzhit) 
      real, intent(in)     :: xyz1(:),xyz2(:)
      logical, intent(out) :: anycross
      real, intent(out)    :: xyzref(:), xyzhit(:)
      end subroutine 

      subroutine d_cart2d_xy(xy, r2)
      real, intent(in)     :: xy(:)
      real, intent(inout)  :: r2(:)
      end subroutine 

      subroutine d_xy2d_cart(xy, r2) 
      real, intent(in)     :: xy(:)
      real, intent(inout)  :: r2(:)
      end subroutine 

      subroutine add_finite_step(xy, dR)
      real, intent(inout) :: xy(:)
      real, intent(in)    :: dR(:)
      end subroutine 
      
      end module
