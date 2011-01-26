ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     NORWECOM grid to sphere conversion
c     ---------------------------------------------------
c     $Revision: 212 $
c     $Date: 2011-01-24 12:08:32 +0100 (Mon, 24 Jan 2011) $
c     $Author: mpay $
c
c     Code taken directly from the NORWECOM model to convert from the NORWECOM grid
c     units to classical lat-lon units
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



C
C***********************************************************************
C
      SUBROUTINE GR2SPH(XGRID, YGRID, XPOLE, YPOLE, DX, DY, YLONG,
     >                  LONG, LATT, IFAIL)
C
C***BEGIN PROLOGUE GR2SPH
C***DATE WRITTEN    7. JAN. 1991
C***REVISION DATE  23. JAN. 1991
C***AUTHOR
C            Bjorn Adlandsvik
C            Institute of Marine Research,
C            Postboks 1870 Nordnes
C            N-5024 Bergen, Norway
C            Email..  bjorn@imr.no
C
C  The code is based on an the routine XYBL from 
C  The Norwegian Meteorological Institute.
C
C***PURPOSE
C  Converts from polar stereographic grid coordinates
C  to spherical coordinates (lattitude and longitude).
C
C***PARAMETERS
C***ON ENTRY
C     XGRID  Real,
C        X grid coordinate.
C     YGRID  Real,
C        Y grid coordinate.
C     XPOLE  Real,
C        X grid coordinate to the North Pole
C     YPOLE  Real,
C        Y grid coordinate to the North Pole.
C     DX     Real,
C        Width in kilometers of grid cell in X direction.
C     DY     Real,
C        Width in kilometers of grid cell in Y direction.
C     YLONG  Real,
C        Longitude in (decimal) degrees of meridian (anti)parallell 
C        to Y-axis.
C
C***ON RETURN
C
C     LONG   Real,
C        Longitude in (decimal) degrees,
C        Normalised to -180 <= LONG <= 180.
C     LATT   Real,
C        Lattitude in (decimal) degrees.
C     IFAIL  Integer,
C        Error indicator, returns
C        2 if DX = 0 or DY = 0,
C        1 if XGRID = XPOLE and YGRID = YPOLE,
C        0 otherwise.
C
C***DESCRIPTION
C  
C  For the coordinate systems see description in routine SPH2GR
C  The North and South Poles does not have well defined longitudes 
C    and are not allowed. This means that 
C    (XGRID, YGRID) = (XPOLE, YPOLE) is not allowed.
C
C  The longitude is computed by
C
C     tan(long-ylong) = (dx * (xgrid-xpole)) / -(dy * (ygrid-ypole))
C
C  and the lattitude by
C
C     tan(45 - 0.5*latt) = r / (rearth * (1 + sin(phinul))
C
C  with
C
C     r**2 = (dx * (xgrid-xp))**2 + (dy * (ygrid-yp))**2
C
C  Normalising is done by finding a rotation number k
C  such that.
C     - 180 <= long + k*360 <= 180
C  that is
C    - long/360 - 1/2 <=  k <= -long/360 + 1/2
C    k is the nearest integer to ( - long/360)
C
C***LOCAL VARIABLES
C
C     LAMBDA		Real,          !MPA: Swapped. This is in fact long
C        Lattitude in radians.    
C     PHI		Real,              !MPA: Swapped. This is in fact lat
C        Longitude in radians.    
C     ALPHA		Real,
C        YLONG in radians.
C     PI 		Real,
C        3.14.....
C     RAD		Real,
C        Conversion factor from degrees to radians.
C     DEG		Real,
C        Conversion factor from radians to degrees.
C     REARTH		Real,
C        Radius (in kilometers) of Earth.
C     PHINUL	   	Real,
C        Lattitude (in radians) of intersecting plane =
C        lattitude of "true" length.
C     R 		Real,
C        Distance (in kilometers) from North Pole in the plane.
C
C***RELATED ROUTINES
C  This routine is the inverse of SPH2GR
C
C***ROUTINES CALLED-NONE
C***END PROLOGUE GR2SPH
C
C  Parameters
C
      REAL XGRID, YGRID, XPOLE, YPOLE, DX, DY, YLONG, LONG, LATT
      INTEGER IFAIL
C
C  Local variables.
C
      REAL LAMBDA, PHI, ALPHA, PI, RAD, DEG, REARTH, PHINUL, R
C
C  Fixed values.  
C
      PARAMETER (PI = 3.14159265)
      PARAMETER (RAD = PI / 180.) 
      PARAMETER (DEG = 1. / RAD)
      PARAMETER (REARTH = 6370.)
      PARAMETER (PHINUL = 60. * RAD)  ! MPA: Commented back in
c.....MPA: Section below commented out
c      INCLUDE 'input.co'
c      IF(ibottom.eq.6)then
c         phinul = 30. * RAD
c      else
c         phinul = 60. * RAD
c      end if


C
C***FIRST EXECUTABLE STATEMENT GR2SPH.
C    
C  Error treatment.
C
      IFAIL = 0
      IF ((XGRID .EQ. XPOLE) .AND. (YGRID .EQ. YPOLE)) THEN
        IFAIL = 1
        WRITE(*,*) '*** ERROR IN GR2SPH -- North Pole is illegal ***'
        GOTO 999
      ENDIF
      IF ((DX .EQ. 0.) .OR. (DY .EQ. 0.)) THEN
        IFAIL = 2
        WRITE(*,*) '*** ERROR IN GR2SPH -- DX or DY is not pos. ***'
        GOTO 999
      ENDIF
C
C  Convert YLONG to radians.
C
      ALPHA = YLONG * RAD
C
C  Longitude.
C 
      LAMBDA = ALPHA + ATAN2(DX*(XGRID-XPOLE), -DY*(YGRID-YPOLE))
C
C  Lattitude.
C
      R = SQRT(DX*DX*(XGRID-XPOLE)*(XGRID-XPOLE) +
     >         DY*DY*(YGRID-YPOLE)*(YGRID-YPOLE))
C
      PHI = .5*PI - 2.*ATAN(R / (REARTH*(1+SIN(PHINUL))))
C
C  Convert to degrees.
C
      LONG = LAMBDA * DEG
      LATT = PHI * DEG
C
C  Normalise longitude.
C
      LONG = LONG - 360. * NINT( LONG/360.)

  999 RETURN
C
C***END GR2SPH
C
      END SUBROUTINE GR2SPH
