ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Physical/mathematical constants
c     ---------------------------------------------------
c     $Rev$
c     $LastChangedDate$
c     $LastChangedBy$ 
c
c     The Earth's equatorial radius =  6378135.0d0 m 
c     The Earth's polar radius       = 6356750.0d0 m
c     The Earth's authalic ("equal area") mean radius                 = 6371005.0761230d0 m 
c     The Earth's volumic radius (radius of a sphere of equal volume) = 6370998.685023d0 m 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module constants
      implicit none
      private                     ! default visibility

      real, parameter,public  :: earth_radius = 6370.0d3                  ! in meters
      real, parameter,public  :: pi           = 3.1415926535897932385d0
      real, parameter,public  :: rad2deg      = 180.0d0/pi                ! factor for conversion from radians to degrees
      real, parameter,public  :: deg2rad      = pi/180.d0 ! factor for conversion from degrees to radians
      real, parameter,public  :: g            = 9.81d0    ! gravity constant
      real, parameter,public  :: rhow0  = 1027.00d00 ! added, changed name rhow -> rhow0, asbjorn
      real, parameter,public  :: gdrho  = g/rhow0 
      end module
