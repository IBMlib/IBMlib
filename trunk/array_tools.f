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
      public ::  histogram_by_bins

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



      subroutine histogram_by_bins(bins, obs, histogram, above, below)
c     ------------------------------------------------------
c     Create a histogram corresponding to dense (adjecent) bins from
c     input data obs.
c     histogram(i) are occurences in obs(:) of numbers
c     in range [bins(i); bins(i+1)]. Entries in obs outside
c     bins are added as occurences in counter (above, below). 
c     The output buffer histogram must have at least (len(bins)-1) entries
c     bins are assumed to be a sorted list of range delimiters
c     ------------------------------------------------------
      real, intent(in)     :: bins(:), obs(:)
      integer, intent(out) :: histogram(:)
      integer, intent(out) :: above, below
      integer              :: n,i,j
c     ------------------------------------------------------
      below     = 0
      above     = 0
      histogram = 0
      n         = size(obs)
      do i = 1,n
         call search_sorted_list(obs(i),bins,j)
         if     (j==0) then
            below = below + 1
         elseif (j==n) then
            above = above + 1
         else
            histogram(j) = histogram(j) + 1
         endif
      enddo
      
c     ------------------------------------------------------
      end subroutine histogram_by_bins


      end module 
