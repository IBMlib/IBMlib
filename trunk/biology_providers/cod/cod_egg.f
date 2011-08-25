      module egg_properties
c     ----------------------------------------------------------
c     Cod egg model 
c     Reference: CJFAS 68 (2011), pp426-443, Table A2
c     ----------------------------------------------------------
      implicit none            
      private  

      public ::  egg_devel_rate
c     ---------------------------------------
                    contains
c     ---------------------------------------


      subroutine egg_devel_rate(temp, rate)
c     ----------------------------------------------------------
c     Assess development progress fraction rate [1/sec]
c     so that expected homeostatic development time is 1/rate
c     at this ambient temperature temp
c     Reference: CJFAS 68 (2011), pp426-443, Table A2
c     ----------------------------------------------------------
      real, intent(in)  :: temp  ! ambient temperature [deg C]
      real, intent(out) :: rate  ! devel rate [1/sec]
      real              :: xpon
c     ----------------------------------------------------------
      xpon = 1.871 - 0.79*log10(temp+2.0)
      rate =  1.0/10**xpon     ! unit = %/day
      rate = rate/100./86400.  ! unit = fraction/sec
      end subroutine egg_devel_rate

      end      ! module
