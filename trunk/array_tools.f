ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Misc generic array tools
c     ---------------------------------------------------
c     $Rev: $
c     $LastChangedDate:  $
c     $LastChangedBy: $  
c   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module array_tools
      implicit none
      private
      public ::  search_sorted_list

c     ---------
      contains
c     ---------      

      subroutine search_sorted_list(v,vlist,i)
c     ------------------------------------------------------
c     Locate the position i in sorted list vlist where v should
c     appear so that  vlist(i) < v < vlist(i+1)
c     Especially return i=0, if  v < vlist(1) and
c     return i=n, if  v > vlist(n)
c     Search by bisection
c     ------------------------------------------------------
      real, intent(in)     :: v, vlist(:)
      integer, intent(out) :: i
      integer              :: ip,itest,n
c     ------------------------------------------------------
      n = size(vlist)
      if (v < vlist(1)) then
         i=0
         return
      endif
      if (v > vlist(n)) then
         i=n
         return
      endif
c     --- do loop: keep vlist(i) < v < vlist(ip) and i<ip
      i  = 1  ! we know v > vlist(1)
      ip = n  ! we know v < vlist(n)
      do while ((ip-i)>1)
         itest = (i+ip)/2 ! integer arithmetics: i<itest<ip
         if (vlist(itest)>v) then
            ip = itest
         else
            i  = itest
         endif
      enddo
c     ------------------------------------------------------
      end subroutine search_sorted_list

      end module 
