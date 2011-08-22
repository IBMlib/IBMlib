      module yolksac_properties
      implicit none            
      private  

      real,parameter :: devel_prefac = 0.27/100./86400. ! daewel etal 2008
      real,parameter :: devel_expo   = 1.495            ! daewel etal 2008

      public ::  yolksac_absorb_rate

c     ---------------------------------------
                    contains
c     ---------------------------------------


      subroutine yolksac_absorb_rate(temp, rate)
c     ----------------------------------------------------------
c     Assess development progress fraction rate [1/sec]
c     so that expected homeostatic development time is 1/rate
c     at this ambient temperature temp
c     ----------------------------------------------------------
      real, intent(in)  :: temp  ! ambient temperature [deg C]
      real, intent(out) :: rate  ! devel rate [1/sec]
c     ----------------------------------------------------------
      rate = devel_prefac * temp**devel_expo
      end subroutine yolksac_absorb_rate

      end      ! module
