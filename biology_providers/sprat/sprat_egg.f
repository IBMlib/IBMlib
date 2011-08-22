      module egg_properties
      implicit none            
      private  

      real,parameter :: devel_prefac = 1.3/100./86400. ! daewel etal 2008
      real,parameter :: devel_expo   = 1.26            ! daewel etal 2008

      public ::  egg_devel_rate
c     ---------------------------------------
                    contains
c     ---------------------------------------


      subroutine egg_devel_rate(temp, rate)
c     ----------------------------------------------------------
c     Assess development progress fraction rate [1/sec]
c     so that expected homeostatic development time is 1/rate
c     at this ambient temperature temp
c     ----------------------------------------------------------
      real, intent(in)  :: temp  ! ambient temperature [deg C]
      real, intent(out) :: rate  ! devel rate [1/sec]
c     ----------------------------------------------------------
      rate = devel_prefac * temp**devel_expo
      end subroutine egg_devel_rate

      end      ! module
