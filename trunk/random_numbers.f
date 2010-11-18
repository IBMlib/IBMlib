ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Random number generators
c     ---------------------------------------------------
c     $Rev$
c     $LastChangedDate$
c     $LastChangedBy$ 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      module random_numbers

      contains

      subroutine get_random_number(r)
c     -------------------------------------------------
c     Proxy function for RNG, until more advanced 
c     generators (like the Merzenne twister) and options
c     for distributions are put in here
c     
c     Generate a random number, 
c     avg(R)=0, Var(R)=1
c     -------------------------------------------------
      real, intent(out) :: r(:)
      real, parameter :: sqrt12 = sqrt(12.0)
      call random_number(r)  ! 0 < r < 1
      r = (r-0.5)*sqrt12
      end subroutine get_random_number
      
      end module
