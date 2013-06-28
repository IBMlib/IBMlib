      module time_services
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This provides time services standardizes time handling between
c     different physical fields. It is a transparent module (public scope)
c     that gives full access to the underlying time_tools module 
c     but hides some internal data structures
c
c     $Rev: $
c     $LastChangedDate:  $
c     $LastChangedBy: $ 
c           
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use time_tools                  ! import clock type
      implicit none
      public           ! default scope (allow fulll access to time_tools)

      public :: get_master_clock
      public :: set_master_clock
      
      interface set_master_clock
        module procedure set_master_clock_i4
        module procedure set_master_clock_from_clock
      end interface
c     -------------------- module data --------------------  
      
      type(clock), target, private :: master_clock
 
c     --------------------
      contains
c     --------------------      


      function   get_master_clock()
c     ------------------------------------------ 
      type(clock), pointer :: get_master_clock
c     ------------------------------------------ 
      get_master_clock => master_clock
c     ------------------------------------------ 
      end function 


      subroutine set_master_clock_i4(idum4)
c     ------------------------------------------
      integer, intent(in) :: idum4(4)
c     ------------------------------------------ 
      call set_clock(master_clock, 
     +               idum4(1), idum4(2), idum4(3), idum4(4))
c     ------------------------------------------ 
      end subroutine 

      subroutine set_master_clock_from_clock(time)
c     ------------------------------------------ 
      type(clock), intent(in) :: time
c     ------------------------------------------ 
      master_clock = time ! local copy
c     ------------------------------------------ 
      end subroutine 

      end module
