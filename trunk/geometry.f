      subroutine cross_2Dline_segments(x0,x1,y0,y1,s,cross)
c     ----------------------------------------------------------
c     Calculate whether the line from x0 to x1 crosses
c     the line from y0 to y1. s is the coordinate along the
c     vector from x0 to x1, and 0<s<1 if cross == .true.
c     If cross == .false. lines do not cross between (x0 to x1)
c     and between (y0 to y1) or lines are parallel.
c
c     Corrected 28 Oct 2010 to ensure only interior solution
c     ----------------------------------------------------------
      implicit none
      real, intent(in)     :: x0(*),x1(*),y0(*),y1(*)
      real, intent(out)    :: s
      logical, intent(out) :: cross
      real                 :: vx(2),vy(2),lvx2,lvy2,det,t
      real,parameter       :: s_parallel = 1.e20
      real,parameter       :: s_undef    = 2.e20
c     ----------------------------------------------------------
      vx    = x1(1:2)-x0(1:2)
      vy    = y1(1:2)-y0(1:2)
      lvx2  = sum(vx*vx)
      lvy2  = sum(vy*vy)
      if ((lvx2 < 1.e-12).or.(lvy2 < 1.e-12)) then ! angle undef
         s = s_undef
         cross=.false.
         return
      endif   
      det = vx(2)*vy(1) - vx(1)*vy(2)     
      if (abs(det)< 1.e-8) then ! lines parallel
         s = s_parallel
         cross=.false.
         return
      endif
c
c     --- now we know lines are not parallel or degenerate ---
c
c     Simplify[Solve[{x01 + vx1 s == y01 + vy1 t, x02 + vx2 s == y02 + vy2 t}, {s,t}]]//FortranForm
c
c              vy2 x01 - vy1 x02 - vy2 y01 + vy1 y02
c         s -> -------------------------------------, 
c                       vx2 vy1 - vx1 vy2
c     
c              vx2 x01 - vx1 x02 - vx2 y01 + vx1 y02
c         t -> -------------------------------------
c                       vx2 vy1 - vx1 vy2
c
c     condition for interior crossing: 0 < s,t < 1
c

      s = (vy(2)*x0(1) - vy(1)*x0(2) - vy(2)*y0(1) + vy(1)*y0(2))/det
      t = (vx(2)*x0(1) - vx(1)*x0(2) - vx(2)*y0(1) + vx(1)*y0(2))/det

      if (in_princip_range(s).and.in_princip_range(t)) then
         cross = .true.
      else
         cross = .false.
      endif
      
      contains

      logical function in_princip_range(x)
      real, intent(in)     :: x
      in_princip_range = ((x >= 0.0).and.(x <= 1.0))
      end function
c     ----------------------------------------------------------
      end subroutine cross_2Dline_segments


c$$$      program jj
c$$$      implicit none
c$$$      real     :: x0(2) = (/0.,  0./)
c$$$      real     :: x1(2) = (/1.,  1./3./)
c$$$      real     :: y0(2) = (/0.,  1./)
c$$$      real     :: y1(2) = (/0.5, 0./) 
c$$$      real     :: s      ! correct answer = 3/7
c$$$      logical  :: cross
c$$$      call cross_2dline_segments(x0,x1,y0,y1,s,cross)
c$$$      write(*,*) cross, s 
c$$$      end program
