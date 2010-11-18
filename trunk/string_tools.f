cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     String tools 
c     ---------------------------------------------------
c     $Rev$
c     $LastChangedDate$
c     $LastChangedBy$ 
c
c     Fortran standard functionality for string handling sucks 
c     
c     Collect a set of usefull string handling utilities
c     to ease life with strings in Fortran
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tokenize(text, start, nwords)
c     -------------------------------------------------
c     Tokenize a text string according to plain spaces
c     (start(i), i=1, nwords) indicates the starting 
c     position in text of each word (block of non-spaces)
c     Return nwords = 0, if no words are resolved    
c     Assume buffer start is sufficently long
c     Unused positions start(nwords+1:) are not defined at exit
c
c     Implementation principle: scan for word starts
c     i.e. the transition " ?", where ? is any non-space
c     character
c     -------------------------------------------------
      character*(*),intent(in) :: text
      integer, intent(out)     :: start(*), nwords
      integer                  :: ltxt, i
c     -------------------------------------------------
      ltxt = len(text)
      nwords = 0
      if (ltxt == 0) return  ! empty text
c
c     ...... handle case if text starts with a non-space ......
c     ...... his is OK, if text only has one character
c
      if (text(1:1) /= " ") then
         nwords = nwords + 1
         start(nwords) = 1
      endif
c
c     ...... for the remaining text look for ......
c     ...... words beginnings " ?"           ......
c
      do i=2,ltxt
         if ((text(i-1:i-1) == " ").and. (text(i:i) /= " ")) then
            nwords = nwords + 1
            start(nwords) = i
         endif
      enddo
      end subroutine tokenize



      subroutine toupper(string)
c     -------------------------------------------------
c     Converts an input string to be all in upper case
c     Useful for string comparisons.
c     Modified from :
c     http://www.star.le.ac.uk/~cgp/fortran.html
c     -------------------------------------------------
      character(len=*), intent(inout) :: string
      integer :: j
c     -------------------------------------------------
      do j = 1,len(string)
        if(string(j:j) >= "a" .and. string(j:j) <= "z") then
            string(j:j) = achar(iachar(string(j:j)) - 32)
        end if
      end do
      end subroutine toupper



      subroutine tolower(string) 
c     -------------------------------------------------
c     Converts an input string to be all in lower case
c     Useful for string comparisons. Modified from toupper() above
c     -------------------------------------------------
      character(len=*), intent(inout) :: string
      integer :: j
c     -------------------------------------------------
      do j = 1,len(string)
        if(string(j:j) >= "A" .and. string(j:j) <= "Z") then
            string(j:j) = achar(iachar(string(j:j)) + 32)
        end if
      end do
      end subroutine tolower
