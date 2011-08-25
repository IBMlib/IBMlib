      module yolksac_properties
c     ----------------------------------------------------------
c     Cod yolk-sac larval model 
c     Reference: CJFAS 68 (2011), pp426-443, Table A2
c     ----------------------------------------------------------

      implicit none            
      private  


      public ::  yolksac_absorb_rate

c     ---------------------------------------
                    contains
c     ---------------------------------------


      subroutine yolksac_absorb_rate(temp, rate)
c     ----------------------------------------------------------
c     Assess development progress fraction rate [1/sec]
c     so that expected homeostatic development time is 1/rate
c     at this ambient temperature temp
c     Reference: CJFAS 68 (2011), pp426-443, Table A2
c     ----------------------------------------------------------
      real, intent(in)  :: temp  ! ambient temperature [deg C]
      real, intent(out) :: rate  ! devel rate [1/sec]
c     ----------------------------------------------------------
      rate = exp(-3.6 + 0.22*temp) ! unit = %/day
      rate = rate/100./86400.      ! unit = fraction/sec

      end subroutine yolksac_absorb_rate

      end      ! module
