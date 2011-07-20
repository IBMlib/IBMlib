ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Limited functionality 1D/2D spline module
c     Based on the plain NR conceptual algorithm of nested splining
c     Extended to provide derivatives through second order as well
c     Flipped nesting order to increase speed compared to original
c    
c     Below, z(:,j) is designated row j and z(i,:) designated column i
c     In fortran the fast index in the first index
c     All subroutines prefixed with int_ are internal to the module
c 
c     Implemented+tested:        ASC June 06, 2011
c     Modularization+overloading:ASC June 21, 2011
c
c     TODO: add optional slopeBC. Default to natural BC
c           For 1D spline this requires yp1,ypn
c           For 2D spline this requires 4 value arrays
c
c     ifort -e90 spline.f
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module spline
      implicit none
      private 

      real,parameter,public :: naturBC     = 1.e30   ! do not change
      real,parameter,public :: naturBC_lim = 0.99e30 ! do not change

      type spline_1D
        real,pointer :: x(:)
        real,pointer :: y(:)
        real,pointer :: y2(:)
      end type
      
      type spline_2D
        real,pointer :: x1(:)   ! size n1
        real,pointer :: x2(:)   ! size n2
        real,pointer :: y(:,:)  ! size n1,n2
        real,pointer :: y2(:,:) ! size n1,n2
      end type

      
      interface setup_spline
         module procedure setup_spline_1D
         module procedure setup_spline_2D
      end interface 

      interface  evaluate_spline
         module procedure evaluate_spline_1D
         module procedure evaluate_spline_2D
      end interface 

      interface delete_spline
         module procedure delete_spline_1D
         module procedure delete_spline_2D
      end interface 
c      
c     ---- define public scope of this module ----
c
      public :: spline_1D
      public :: spline_2D
      public :: setup_spline
      public :: evaluate_spline
      public :: delete_spline

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                       contains
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      



      SUBROUTINE setup_spline_1D(spli1D,x,y)
c     ---------------------------------------
      type(spline_1D), intent(out) :: spli1D
      REAL,intent(in)              :: x(:),y(:)
      integer                      :: n
      REAL                         :: yp1,ypn
c     ---------------------------------------
      yp1 = naturBC
      ypn = naturBC
      n = size(x) ! assert = size(y) 
      allocate(spli1D%x(n))
      allocate(spli1D%y(n))
      allocate(spli1D%y2(n))
      spli1D%x = x ! copy data
      spli1D%y = y ! copy data
      call int_setup_spline(spli1D%x,spli1D%y,yp1,ypn,spli1D%y2)
c     ---------------------------------------
      end SUBROUTINE setup_spline_1D



      SUBROUTINE setup_spline_2D(spli2D,x1,x2,y)
c     ---------------------------------------
      type(spline_2D), intent(out) :: spli2D
      REAL,intent(in)              :: x1(:),x2(:),y(:,:)
      integer                      :: n1,n2
c     ---------------------------------------
      n1 = size(x1) ! assert = size(y,1) 
      n2 = size(x2) ! assert = size(y,2) 
      allocate(spli2D%x1(n1))
      allocate(spli2D%x2(n2))
      allocate(spli2D%y(n1,n2))
      allocate(spli2D%y2(n1,n2))
      spli2D%x1 = x1 ! copy data
      spli2D%x2 = x2 ! copy data
      spli2D%y  = y  ! copy data
      call int_setup_2D_spline(spli2D%x1,spli2D%x2,spli2D%y,spli2D%y2)
c     ---------------------------------------
      end SUBROUTINE setup_spline_2D



      SUBROUTINE delete_spline_1D(spli1D)
c     ---------------------------------------
      type(spline_1D), intent(out) :: spli1D
c     ---------------------------------------
      deallocate(spli1D%x)
      deallocate(spli1D%y)
      deallocate(spli1D%y2)
c     ---------------------------------------
      end SUBROUTINE delete_spline_1D


      SUBROUTINE delete_spline_2D(spli2D)
c     ---------------------------------------
      type(spline_2D), intent(out) :: spli2D
c     ---------------------------------------
      deallocate(spli2D%x1)
      deallocate(spli2D%x2)
      deallocate(spli2D%y)
      deallocate(spli2D%y2)
c     ---------------------------------------
      end SUBROUTINE delete_spline_2D


      
      SUBROUTINE evaluate_spline_1D(spli1D,x,deriv,y)
c     ---------------------------------------
      type(spline_1D), intent(in) :: spli1D
      REAL,intent(in)             :: x
      integer,intent(in)          :: deriv
      REAL,intent(out)            :: y
c     ---------------------------------------
      call int_evaluate_spline(spli1D%x,spli1D%y,spli1D%y2,x,deriv,y)
c     ---------------------------------------
      end SUBROUTINE evaluate_spline_1D



      SUBROUTINE evaluate_spline_2D(spli2D,x1,x2,deriv,y)
c     ---------------------------------------
      type(spline_2D), intent(in) :: spli2D
      REAL,intent(in)             :: x1,x2
      integer,intent(in)          :: deriv(2)
      REAL,intent(out)            :: y
c     ---------------------------------------
      call int_evaluate_2D_spline(spli2D%x1,spli2D%x2,
     +                         spli2D%y,spli2D%y2,
     +                         x1,x2,deriv,y)
c     ---------------------------------------
      end SUBROUTINE evaluate_spline_2D


cccccccccccccccccccccccc internal subroutines cccccccccccccccccccccccc



      SUBROUTINE int_setup_2D_spline(x1a,x2a,ya,y2a)
c     ----------------------------------------------- 
c     Setup 2d spline as rowwise splines
c   
c     Input data is ya(:,:) on a regular grid and
c     x1a(:) corresponds to grid along first index and
c     x2a(:) corresponds to grid along second index
c     Intermediate output is the matrix of second derivatives y2a 
c     In principle x2a is not needed to be included 
c     (because splining is row wise), but included it to 
c     allow nesting flip or so later, if needed
c     -----------------------------------------------
      REAL,intent(in)   :: x1a(:),x2a(:),ya(:,:)
      REAL,intent(out)  :: y2a(:,:)
      INTEGER           :: k
      do k=1,size(ya,2)
         call int_setup_spline(x1a,ya(:,k),naturBC,naturBC,y2a(:,k))
      enddo
c
      END SUBROUTINE int_setup_2D_spline



      SUBROUTINE int_evaluate_2D_spline(x1a,x2a,ya,y2a,x1,x2,deriv,y)
c     -----------------------------------------------
c     Evaluate 2d spline in point (x1,x2) in scalar y 
c
c     Input data is ya(:,:) on a regular grid and
c     x1a(:) corresponds to grid along first index and
c     x2a(:) corresponds to grid along second index
c     and y2a is the matrix of second derivatives 
c     generated by setup_2D_spline 
c     deriv(i,j) refers to derivative order along 
c     first/second index respectively. deriv(0,0) gives
c     value spline in point (x1,x2)
c     -----------------------------------------------
      REAL,intent(in)    :: x1a(:),x2a(:),ya(:,:),y2a(:,:),x1,x2
      INTEGER,intent(in) :: deriv(2)
      REAL,intent(out)   :: y
      INTEGER            :: k
      REAL               :: ytmp(size(ya,2))  ! automatic var
      REAL               :: y2tmp(size(ya,2)) ! automatic var
c     ---------------------------------
c     create tmp vertial array to splined values ytmp(:) at x1 
c     corresponding to grid x2a(:)     
      do k=1,size(ya,2)
         call int_evaluate_spline(x1a,ya(:,k),y2a(:,k),x1,
     +                            deriv(1),ytmp(k))
      enddo
c     --- setup intermediate vertical spline ---
      call int_setup_spline(x2a,ytmp,naturBC,naturBC,y2tmp)
c     --- evaluate intermediate vertical spline at x2 ---
      call int_evaluate_spline(x2a,ytmp,y2tmp,x2,deriv(2),y)
c
      END SUBROUTINE int_evaluate_2D_spline




      SUBROUTINE int_setup_spline(x,y,yp1,ypn,y2)
c     -----------------------------------------------
c     Setup 1d spline    
c
c     Input data is ya(:) on a grid xa(:) and 
c     yp1,ypn are slopes of y in  xa(1) and xa(last) resplectively
c     Output is the vector y2a(:) of second derivatives on
c     same grid xa(:) as needed by splint
c     Based on common implementation spline
c     -----------------------------------------------
      REAL,intent(in)     :: x(:),y(:),yp1,ypn
      REAL,intent(out)    :: y2(:)
      INTEGER i,k,n
      REAL p,qn,sig,un,u(size(x))
c     -----------------------------------------------
      n = size(x)
      if (yp1.gt.naturBC_lim) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
      enddo

      if (ypn.gt.naturBC_lim) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      enddo
      
      END SUBROUTINE int_setup_spline



      SUBROUTINE int_evaluate_spline(xa,ya,y2a,x,deriv,y)
c     -----------------------------------------------
c     Evaluate 1d spline
c  
c     Input data is ya(:) on a grid xa(:) and 
c     vector y2a(:) of second derivatives as generated by spline
c     deriv = 0 -> return splined value in y
c     deriv = 1 -> return splined slope in y
c     deriv = 2 -> return splined second derivatives in y
c     Notice that second derivative of a cubic spline is
c     continuous everywhere and piecewise linear.
c
c     Based on common implementation splint
c     -----------------------------------------------
      REAL,intent(in)     :: xa(:),ya(:),y2a(:),x
      INTEGER,intent(in)  :: deriv
      REAL,intent(out)    :: y

      INTEGER n,k,khi,klo
      REAL a,b,h
      n  = size(xa)
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) stop 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      if (deriv==0) then              ! see Numerical Recipes p.114
         y = a*ya(klo) + b*ya(khi) + ((a**3-a)*y2a(klo) 
     +               + (b**3-b)*y2a(khi))*(h**2)/6.
      elseif (deriv==1) then          ! see Numerical Recipes p.114
         y = (ya(khi)-ya(klo))/h 
     +       - (3.*a**2 - 1.)*h*y2a(klo)/6.
     +       + (3.*b**2 - 1.)*h*y2a(khi)/6.
      elseif (deriv==2) then          ! see Numerical Recipes p.114
         y = a*y2a(klo) + b*y2a(khi)  
      else
         stop "unknown deriv argument"
      endif
      
      END SUBROUTINE int_evaluate_spline

      end ! module


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccc         T E S T     B L O C K   ccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$
c$$$      program test
c$$$      use spline
c$$$      implicit none 
c$$$
c$$$      integer, parameter :: n1 = 21, n2 =31
c$$$      real               :: x1a(n1),x2a(n2)
c$$$      real               :: ya(n1,n2),y2a(n1,n2)
c$$$      real               :: y2n1(n1),yan1(n1),x,y,dydx,dy2dx2
c$$$      real               :: y00, ym0, yp0, y0m, y0p, ymm, ymp, ypm, ypp
c$$$      real               :: d10_y, d01_y, d11_y, d20_y, d02_y
c$$$      real               :: dx1, dx2, x1, x2
c$$$      integer            :: i, j, d(2)
c$$$c
c$$$c     ------------ setup data ------------
c$$$c      
c$$$      do i=1,n1
c$$$         x1a(i) = i-1
c$$$      enddo
c$$$      x1a = 2 + 3*x1a   ! [2; 3*(n1-1)=60]
c$$$      do j=1,n2
c$$$         x2a(j) = j-1
c$$$      enddo
c$$$      x2a = 1 + 0.5*x2a ! [1; 0.5*(n2-1)=16]
c$$$      do i=1,n1
c$$$         x = x1a(i)
c$$$         do j=1,n2
c$$$            y = x2a(j)
c$$$            ya(i,j) = sin(x) + cos(1.441*x)*sqrt(y)
c$$$         enddo
c$$$      enddo
c$$$c
c$$$c     ------------ test 1D spline ------------
c$$$c     
c$$$      yan1 = ya(:,4)
c$$$      call setup_spline(x1a,ya(:,4),naturBC,naturBC,y2n1)
c$$$      do x=x1a(1),x1a(n1),0.01
c$$$         call evaluate_spline(x1a,yan1,y2n1,x,0,y)
c$$$         call evaluate_spline(x1a,yan1,y2n1,x,1,dydx)
c$$$         call evaluate_spline(x1a,yan1,y2n1,x,2,dy2dx2)
c$$$         write(12,*) x,y,dydx,dy2dx2
c$$$      enddo
c$$$c
c$$$c     ------------ test 2D spline ------------
c$$$c       
c$$$      call setup_2D_spline(x1a,x2a,ya,y2a)
c$$$      x1 = x1a(4) + 0.37
c$$$      x2 = x2a(9) + 0.182
c$$$c     --- splined derivatives ---      
c$$$      d = (/1, 0/)
c$$$      call evaluate_2D_spline(x1a,x2a,ya,y2a,x1,x2,d,d10_y)
c$$$      d = (/0, 1/)
c$$$      call evaluate_2D_spline(x1a,x2a,ya,y2a,x1,x2,d,d01_y)
c$$$      d = (/1, 1/)
c$$$      call evaluate_2D_spline(x1a,x2a,ya,y2a,x1,x2,d,d11_y)
c$$$      d = (/0, 2/)
c$$$      call evaluate_2D_spline(x1a,x2a,ya,y2a,x1,x2,d,d02_y)
c$$$      d = (/2, 0/)
c$$$      call evaluate_2D_spline(x1a,x2a,ya,y2a,x1,x2,d,d20_y)
c$$$c     --- numerical derivatives --- 
c$$$      d   = (/0, 0/)
c$$$      dx1 = 0.05
c$$$      dx2 = 0.05
c$$$      call evaluate_2D_spline(x1a,x2a,ya,y2a,x1    ,x2    ,d,y00)
c$$$      call evaluate_2D_spline(x1a,x2a,ya,y2a,x1-dx1,x2    ,d,ym0)
c$$$      call evaluate_2D_spline(x1a,x2a,ya,y2a,x1+dx1,x2    ,d,yp0)
c$$$      call evaluate_2D_spline(x1a,x2a,ya,y2a,x1    ,x2-dx2,d,y0m)
c$$$      call evaluate_2D_spline(x1a,x2a,ya,y2a,x1    ,x2+dx2,d,y0p)
c$$$      call evaluate_2D_spline(x1a,x2a,ya,y2a,x1-dx1,x2-dx2,d,ymm)
c$$$      call evaluate_2D_spline(x1a,x2a,ya,y2a,x1-dx1,x2+dx2,d,ymp)
c$$$      call evaluate_2D_spline(x1a,x2a,ya,y2a,x1+dx1,x2-dx2,d,ypm)
c$$$      call evaluate_2D_spline(x1a,x2a,ya,y2a,x1+dx1,x2+dx2,d,ypp)
c$$$c
c$$$      write(*,333) "dy_dx1    ", d10_y, (yp0-ym0)/2/dx1
c$$$      write(*,333) "dy_dx2    ", d01_y, (y0p-y0m)/2/dx2
c$$$      write(*,333) "dy2_dx1dx2", d11_y, (ypp+ymm-ymp-ypm)/2/dx1/2/dx2
c$$$      write(*,333) "dy2_dx1dx1", d20_y, (yp0+ym0-2*y00)/dx1/dx1
c$$$      write(*,333) "dy2_dx2dx2", d02_y, (y0p+y0m-2*y00)/dx2/dx2
c$$$c
c$$$ 333  format(a," : splin=", f12.7, " ndiff=", f12.7)
c$$$      end program

