c     ========= TODO extend with "deriv" argument
c                    0(value), 3(zgrad), 12(horiz grad), 123(full 3Dgrad)
c                    upgrade to module and have an overloaded interface, 
c                    accepting scalar/vector output request
c                    currently, below, only scalar output enabled


      subroutine interp_2Dbox_data(sx,sy,vc,deriv,result)
c     ------------------------------------------------------------
c     Basic bilinear interpolation of 2D grid data
c     deriv = 0 gives value, deriv = (1,2) gives derivative (wrt. sx,sy) along (x,y)
c     Currently, only deriv = 0 is implemented, until other derivatives 
c     are needed.
c     ------------------------------------------------------------
      implicit none
      real,intent(in)  :: sx,sy
      real,intent(in)  :: vc(*) ! 4 corners values : (v00,v01,v10,v11)
      integer,intent(in) :: deriv
      real,intent(out) :: result
      real             :: v00,v01,v10,v11, vl,vu
c     ------------------------------------------------------------
      v00 = vc(1)
      v01 = vc(2)
      v10 = vc(3)
      v11 = vc(4)

      vl  = (1-sx)*v00 + sx*v10
      vu  = (1-sx)*v01 + sx*v11

      if     (deriv == 0) then
         result = (1-sy)*vl + sy*vu
      else
         stop "interp_2Dbox_data: unhandled deriv request"
      endif
c     ------------------------------------------------------------
      end subroutine interp_2Dbox_data


      subroutine interp_irregbox_data(sx,sy,z,zc,vc,deriv,result,istat)
c     ------------------------------------------------------------
c     Interpolate between 8 corner (zc) value points (vc) of an 
c     irregular box at coordinate (sx,sy,z), where all pillars 
c     (vertical sides) are strictly vertical 
c     (so vertical projection is a rhomb)
c     zc is arranged as (z000, z001, ... z111) and similar for vc
c     where the integer indices refres to xyz directions
c 
c     Use generalized trilinear interpolation
c     (sx,sy) are rhomb coordinates on vertical projection,
c     0 < sx,sy < 1, and z is the vertical coordinate
c
c     deriv = 0 gives value, deriv = (1,2,3) gives derivative (wrt. sx,sy,z) along (x,y,z)
c     Currently, only deriv = (0,3) is implemented, until other derivatives 
c     are needed.
c
c     istat = 0: regular interpolation performed
c     istat = 1: singular vertical - performed mid point interpolation
c     ------------------------------------------------------------
      implicit none
      real,intent(in)     :: sx,sy,z
      real,intent(in)     :: zc(*) ! 8 corners:         (z000,z001, ... z111)
      real,intent(in)     :: vc(*) ! 8 corners values : (v000,v001, ... v111)
      integer,intent(in)  :: deriv
      real,intent(out)    :: result
      integer,intent(out) :: istat

      real             :: zl,zu,vl,vu,sz,dz
      real             :: z000,z001,z010,z011,z100,z101,z110,z111
      real             :: v000,v001,v010,v011,v100,v101,v110,v111
      real, parameter  :: very_small = 1.e-12
c     ------------------------------------------------------------
      z000 = zc(1)
      z001 = zc(2)
      z010 = zc(3)
      z011 = zc(4)
      z100 = zc(5)
      z101 = zc(6)
      z110 = zc(7)
      z111 = zc(8)
c
      v000 = vc(1)
      v001 = vc(2)
      v010 = vc(3)
      v011 = vc(4)
      v100 = vc(5)
      v101 = vc(6)
      v110 = vc(7)
      v111 = vc(8)
c     ------- upper/lower lid horizontal interpolation at (sx,sy) 
c             (zl,vl) corresponds to sz=0
c             (zu,vl) corresponds to sz=1  (zu>zl))
      zl = (1.0-sy)*((1-sx)*z000 + sx*z100) +
     +          sy *((1-sx)*z010 + sx*z110) 
      zu = (1.0-sy)*((1-sx)*z001 + sx*z101) +
     +          sy *((1-sx)*z011 + sx*z111) 
c
      vl = (1.0-sy)*((1-sx)*v000 + sx*v100) +
     +          sy *((1-sx)*v010 + sx*v110) 
      vu = (1.0-sy)*((1-sx)*v001 + sx*v101) +
     +          sy *((1-sx)*v011 + sx*v111) 
c     ------- vertical interpolation at z 
      dz = zu-zl
      if (abs(dz) < very_small) then
         sz    = 0.5
         istat = 1
      else ! normal case
         sz = (z-zl)/(zu-zl)
         istat = 0
      endif

c$$$      if ((sz<0.).or.(sz>1.)) then
c$$$         write(*,*) "zc=",zc(1:8)  
c$$$         write(*,*) "vc=",vc(1:8)
c$$$         write(*,*) "zl,zu=",zl,zu
c$$$         write(*,*) "z, sz=",z, sz
c$$$         stop "interpolate_irregbox_data: illegal coordinate"
c$$$      endif

      if     (deriv == 0) then
         result = vl + sz*(vu-vl)
      elseif (deriv == 3) then
         if (istat == 0) then
            result = (vu-vl)/dz ! normal case
         else
            result = 0  ! zero derivative in this case
         endif   
      else
         stop "interp_irregbox_data: unhandled deriv request"
      endif

      end subroutine interp_irregbox_data



