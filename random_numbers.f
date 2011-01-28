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
 
      interface get_random_number
        module procedure get_random_number_vec
        module procedure get_random_number_single
      end interface

      contains

      subroutine get_random_number_vec(r)
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
      end subroutine get_random_number_vec

      
      subroutine get_random_number_single(r)
c     -------------------------------------------------
c     Overloaded version that returns just a single numbe
c     rather than a vector
c     -------------------------------------------------
      real, intent(out) :: r
      real :: r1(1)
      call get_random_number_vec(r1)  
      r=r1(1)
      end subroutine get_random_number_single

      end module
