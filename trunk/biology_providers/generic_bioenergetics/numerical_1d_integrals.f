      module numerical_1d_integrals
      implicit none   
      private

      interface trapez_integral
        module procedure trapez_integral0
        module procedure trapez_integral1
        module procedure trapez_integral2
      end interface
      public :: trapez_integral

      contains

      subroutine trapez_integral0(y, h, intg) 
c     ------------------------------------------------------------------
c     Evaluate the trapez integral intg on same regular sampling grid. 
c     ------------------------------------------------------------------
      real,intent(in)  :: y(:)
      real,intent(in)  :: h
      real,intent(out) :: intg 
      integer          :: n,i
      real             :: w(size(y))
      real             :: z
c     ------------------------------------------------------------------
      n      = size(y)
      w      = 1.0
      w(1)   = 0.5
      w(n)   = 0.5   
c
      z         = 0.0      
      do i=1,n
         z         = z + w(i)*y(i)          
      enddo
      intg         = h*z    
      end subroutine trapez_integral0


      subroutine trapez_integral1(y, dy, h, intg, dintg_dx0, dintg_dx1) 
c     ------------------------------------------------------------------
c     Evaluate the trapez integral intg and the first order end points derivatives
c     (dintg_dx0, dintg_dx1) of the integrand y(x) on a regular sampling grid
c     x = [x0, x0+h, ..., x1](not in argument list - not needed explicitly) 
c     with regular grid spacing h. dy is the derivative of y, dydx(x), 
c     on same regular sampling grid.
c     
c     The end point derivatives (dintg_dx0, dintg_dx1) are consistent with
c     values obtained by numerical differentiation when using a 
c     finite sampling grid x. When a finite sampling grid is used, this
c     is slightly different from the analytic end point rule. When 
c     sampling grid x becomes infinitely dense, (dintg_dx0, dintg_dx1) 
c     should converge to the result of the analytic end point rule
c
c     Validated by numerical differentiation    
c     ------------------------------------------------------------------
      real,intent(in)  :: y(:)
      real,intent(in)  :: dy(:)
      real,intent(in)  :: h
      real,intent(out) :: intg                  ! trapez integral
      real,intent(out) :: dintg_dx0, dintg_dx1  ! 1st derivatives
      integer          :: n,i
      real             :: w(size(y))
      real             :: z,dz_dx0,dz_dx1 
      real             :: dh_dx0, dh_dx1, ds_dx0, ds_dx1
c     ------------------------------------------------------------------
      n      = size(y)
      w      = 1.0
      w(1)   = 0.5
      w(n)   = 0.5
      dh_dx0 = -1.0/(n-1)  ! higher derivatives of h is zero
      dh_dx1 =  1.0/(n-1)  ! higher derivatives of h is zero
c
      z         = 0.0 
      dz_dx0    = 0.0
      dz_dx1    = 0.0
      do i=1,n
         ds_dx0    = 1.0*(n-i)/(n-1) ! higher derivatives of s is zero
         ds_dx1    = 1.0*(i-1)/(n-1) ! higher derivatives of s is zero
         z         = z + w(i)*y(i)
         dz_dx0    = dz_dx0 + w(i)*dy(i)*ds_dx0
         dz_dx1    = dz_dx1 + w(i)*dy(i)*ds_dx1  
      enddo
      intg         = h*z
      dintg_dx0    = h*dz_dx0 + dh_dx0*z
      dintg_dx1    = h*dz_dx1 + dh_dx1*z           
      end subroutine trapez_integral1



      subroutine trapez_integral2(y, dy, d2y, h, intg, 
     +                       dintg_dx0, dintg_dx1,
     +                       d2intg_dx0x0, d2intg_dx1x1, d2intg_dx0x1) 
c     ------------------------------------------------------------------
c     Evaluate the trapez integral intg and the first order  (dintg_dx0, dintg_dx1)
c     and second order end points derivatives (d2intg_dx0x0, d2intg_dx1x1, d2intg_dx0x1)
c     of the integrand y(x) on a regular sampling grid
c     x = [x0, x0+h, ..., x1](not in argument list - not needed explicitly) 
c     with regular grid spacing h. dy is the derivative of y, dydx(x), 
c     and d2y the second derivative of y, d2ydx2(x), all on same regular sampling grid.
c     
c     The end point derivatives (dintg_dx0, dintg_dx1) and
c     (d2intg_dx0x0, d2intg_dx1x1, d2intg_dx0x1) are consistent with
c     values obtained by numerical differentiation when using a 
c     finite sampling grid x. When a finite sampling grid is used, this
c     is slightly different from the analytic end point rule. When 
c     sampling grid x becomes infinitely dense (dintg_dx0, dintg_dx1) and
c     (d2intg_dx0x0, d2intg_dx1x1, d2intg_dx0x1) 
c     should converge to the result of the analytic end point rule
c
c     Validated by numerical differentiation    
c     ------------------------------------------------------------------
      real,intent(in)  :: y(:)
      real,intent(in)  :: dy(:)
      real,intent(in)  :: d2y(:)
      real,intent(in)  :: h
      real,intent(out) :: intg                  ! trapez integral
      real,intent(out) :: dintg_dx0, dintg_dx1  ! 1st derivatives
      real,intent(out) :: d2intg_dx0x0, d2intg_dx1x1, d2intg_dx0x1 ! 2nd derivatives
      integer          :: n,i
      real             :: w(size(y))
      real             :: z,dz_dx0,dz_dx1,d2z_dx0x0,d2z_dx0x1,d2z_dx1x1
      real             :: dh_dx0, dh_dx1, ds_dx0, ds_dx1
c     ------------------------------------------------------------------
      n      = size(y)
      w      = 1.0
      w(1)   = 0.5
      w(n)   = 0.5
      dh_dx0 = -1.0/(n-1)  ! higher derivatives of h is zero
      dh_dx1 =  1.0/(n-1)  ! higher derivatives of h is zero
c
      z         = 0.0 
      dz_dx0    = 0.0
      dz_dx1    = 0.0
      d2z_dx0x0 = 0.0
      d2z_dx1x1 = 0.0
      d2z_dx0x1 = 0.0
      do i=1,n
         ds_dx0    = 1.0*(n-i)/(n-1) ! higher derivatives of s is zero
         ds_dx1    = 1.0*(i-1)/(n-1) ! higher derivatives of s is zero
         z         = z + w(i)*y(i)
         dz_dx0    = dz_dx0 + w(i)*dy(i)*ds_dx0
         dz_dx1    = dz_dx1 + w(i)*dy(i)*ds_dx1
         d2z_dx0x0 = d2z_dx0x0 + w(i)*d2y(i)*ds_dx0*ds_dx0
         d2z_dx1x1 = d2z_dx1x1 + w(i)*d2y(i)*ds_dx1*ds_dx1
         d2z_dx0x1 = d2z_dx0x1 + w(i)*d2y(i)*ds_dx0*ds_dx1
      enddo
      intg         = h*z
      dintg_dx0    = h*dz_dx0 + dh_dx0*z
      dintg_dx1    = h*dz_dx1 + dh_dx1*z
      d2intg_dx0x0 = h*d2z_dx0x0 + 2.0*dh_dx0*dz_dx0 
      d2intg_dx1x1 = h*d2z_dx1x1 + 2.0*dh_dx1*dz_dx1 
      d2intg_dx0x1 = h*d2z_dx0x1 + dh_dx0*dz_dx1 + dh_dx1*dz_dx0 
      end subroutine trapez_integral2

      end module


c      program test
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc     ifort -e90 test_trapez_integrals.f; a.out
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      use numerical_1d_integrals
c      implicit none   
c      
c      integer, parameter :: n=30
c      real               :: y(n),dy(n),d2y(n),x0,x1,dx,xi
c      integer            :: i
c      real               :: t00,h,tm0,tp0,t0m,t0p,tmm,tmp,tpm,tpp,t
c      real               :: nd01,nd10,nd02,nd20,nd11
c      real               :: dt_dx0,dt_dx1,d2t_dx0x0,d2t_dx1x1,d2t_dx0x1
c      x0 = 1.13
c      x1 = 5.28
c      dx = 0.01
c      
c      call set_arrays(x0-dx,x1   ,h,y,dy,d2y)
c      call trapez_integral(y, h, tm0) 
c      call set_arrays(x0+dx,x1   ,h,y,dy,d2y)
c      call trapez_integral(y, h, tp0)
c      call set_arrays(x0   ,x1-dx,h,y,dy,d2y)
c      call trapez_integral(y, h, t0m) 
c      call set_arrays(x0   ,x1+dx,h,y,dy,d2y)
c      call trapez_integral(y, h, t0p) 
cc
c      call set_arrays(x0-dx,x1-dx,h,y,dy,d2y)
c      call trapez_integral(y, h, tmm) 
c      call set_arrays(x0-dx,x1+dx,h,y,dy,d2y)
c      call trapez_integral(y, h, tmp) 
c      call set_arrays(x0+dx,x1-dx,h,y,dy,d2y)
c      call trapez_integral(y, h, tpm) 
c      call set_arrays(x0+dx,x1+dx,h,y,dy,d2y)
c      call trapez_integral(y, h, tpp) 
c      call set_arrays(x0   ,x1   ,h,y,dy,d2y)
c      call trapez_integral(y, h, t00)           ! base as final case
cc
c      nd01 = (t0p-t0m)/2./dx
c      nd10 = (tp0-tm0)/2./dx
c      nd02 = (t0p+t0m-2*t00)/dx**2
c      nd20 = (tp0+tm0-2*t00)/dx**2
c      nd11 = (tpp+tmm-tmp-tpm)/4./dx**2
cc
c      write(*,*) "test trapez_integral0: t00=",t00
c      call trapez_integral(y, dy, h, t,dt_dx0,dt_dx1)
c      write(*,*) "test trapez_integral1: t     =",t
c      write(*,*) "test trapez_integral1: dt_dx0=",dt_dx0,nd10
c      write(*,*) "test trapez_integral1: dt_dx1=",dt_dx1,nd01
c      call trapez_integral(y, dy, d2y, h, t,dt_dx0,dt_dx1,
c     +                     d2t_dx0x0, d2t_dx1x1,d2t_dx0x1)
c      write(*,*) "test trapez_integral2: t     =",t
c      write(*,*) "test trapez_integral2: dt_dx0=",dt_dx0,nd10
c      write(*,*) "test trapez_integral2: dt_dx1=",dt_dx1,nd01
c      write(*,*) "test trapez_integral2: d2t_dx0x0=",d2t_dx0x0,nd20
c      write(*,*) "test trapez_integral2: d2t_dx1x1=",d2t_dx1x1,nd02
c      write(*,*) "test trapez_integral2: d2t_dx0x1=",d2t_dx0x1,nd11
c
c      contains
c      
c      subroutine set_arrays(x0,x1,h,y,dy,d2y)
c      real,intent(in)  :: x0,x1
c      real,intent(out) :: h,y(n),dy(n),d2y(n)
c      real             :: xi
c      integer          :: i
c      h  = (x1-x0)/(n-1.0)
c      do i=1,n
c         xi = x0 + (x1-x0)*(i-1.0)/(n-1.0)
c         y(i)   = sin(xi**2)
c         dy(i)  = 2*xi*cos(xi**2)
c         d2y(i) = 2*cos(xi**2) - 4*xi*xi*sin(xi**2)
c      enddo
c      end subroutine set_arrays
c
c      end
