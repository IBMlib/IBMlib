ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     NORWECOM sphere to grid conversion
c     ---------------------------------------------------
c     $Revision: 212 $
c     $Date: 2011-01-24 12:08:32 +0100 (Mon, 24 Jan 2011) $
c     $Author: mpay $
c
c     Code taken directly from the NORWECOM model to convert from spherical (lat-lon)
c     units to the NORWECOM grid
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C***********************************************************************
C
      SUBROUTINE SPH2GR(LONG, LATT, XPOLE, YPOLE, DX, DY, YLONG,
     >                  XGRID, YGRID, IFAIL)
C
C***BEGIN PROLOGUE SPH2GR
C***DATE WRITTEN   910104 (YYMMDD)
C***REVISION DATE  910104 (YYMMDD)
C***AUTHOR
C            Bjorn Adlandsvik
C            Institute of Marine Research,
C            Postboks 1870 Nordnes
C            N-5024 Bergen, Norway
C            Email..  bjorn@imr.no
C
C  The code is based on an the routine BLXY from 
C  The Norwegian Meteorological Institute.
C
C***PURPOSE
C  Converts from spherical coordinates (lattitude and longitude)
C  to polar stereographic grid coordinates.
C
C***PARAMETERS
C***ON ENTRY
C     LONG   Real,	
C  	 Longitude in (decimal) degrees.
C     LATT   Real, 
C        Lattitude in (decimal) degrees,
C        Requirement: -90. < LATT <= 90.    
C     XPOLE  Real,
C        X grid coordinate of the North Pole
C     YPOLE  Real,
C        Y grid coordinate of the North Pole.
C     DX     Real,
C        Width in kilometers of grid cell in X direction.
C     DY     Real,
C        Width in kilometers of grid cell in Y direction.
C     YLONG  Real,
C        Longitude in (decimal) degrees of meridian (anti)parallell 
C        to the Y-axis. 
C
C***ON RETURN
C
C     XGRID  Real,
C        X grid coordinate. 
C     YGRID  Real,
C        Y grid coordinate.
C     IFAIL  Integer,
C        Error indicator, returns
C        2 if DX <= 0 or DY <= 0,
C        1 if illegal value of LATT,
C        0 otherwise.
C
C***DESCRIPTION
C  
C  Spherical coordinates:
C  These are ordinary longitude and lattitude in decimal degrees. 
C  Decimal degrees means that the angle is described as a decimal 
C     number, i.e. 64.5 is used for 64 degrees and 30 minutes.
C  Positive values of longitude means easterly and negative westerly.
C     All values are allowed and are treated modulo 360.
C  Positive values of lattitude means northern and negative southern.
C     The values are restricted to the intervall (-90., 90].
C     The South Pole is not included since the polar stereographic
C     projection is not defined there. The North Pole is treated
C     correctly with any value of longitude.
C 
C  Grid coordinates:
C  The grid coordinate system is based on the polar stereographic
C     projection, "true" at 60 degrees north. This longitude can be 
C     changed by varying the constant PHINUL below.
C  The grids origin in the plane, its orientation and the scaling
C     of the axes are given by XP,YP,YLONG,DX,DY above.
C  The grid coordinates of the center of the (i, j)-th grid cell
C     are (i-0.5, j-0.5)
C
C  The projected radius r is computed by the formula
C     r = rearth * cos(phi) * (1 + sin(phinul)) / (1 + sin(phi))
C
C  Then the grid coordinate are computed as follows
C     xgrid = xpole + r * sin(lambda-alpha) / dx
C     ygrid = ypole - r * cos(lambda-alpha) / dy
C
C***LOCAL VARIABLES
C
C     LAMBDA		Real        !MPA: Swapped. This is in fact lon
C        lattitude in radians.
C     PHI		Real            !MPA: Swapped. This is in fact lat
C        Longitude in radians.
C     ALPHA		Real
C        YLONG in radians.
C     RAD		Real
C        Conversion factor from degrees to radians.
C     REARTH		Real
C        Radius (in kilometers) of Earth.
C     PHINUL
C        Lattitude (in radians) of intersecting plane =
C        lattitude of "true" length.
C     R 		Real
C        Distance (in kilometers) from North Pole in the plane.
C
C***REMARKS
C  Parameters from some often used grids : 
C                           xpole, ypole,  dx,  dy, ylong
C  15570 grid :              181., 128.5, 20., 20.,  58.
C  Barents Sea grid (95,95):  90., 128.5, 20., 20.,  58
C  Hindcast grid :            11.,  65. , 75., 75., -32.
C  North Sea grid            352., 256. , 10., 10.,  58.
C
C***RELATED ROUTINES
C  This routine is the inverse of GR2SPH
C
C***ROUTINES CALLED-NONE
C***END PROLOGUE SPH2GR
C
C  Parameters
C
      REAL LONG, LATT, XPOLE, YPOLE, DX, DY, YLONG, XGRID, YGRID
      INTEGER IFAIL
C
C  Local variables.
C
      REAL LAMBDA, PHI, ALPHA, RAD, REARTH, PHINUL, R 
C
C  Fixed values.  
C
      PARAMETER (RAD = 3.14159265 / 180.)
      PARAMETER (REARTH = 6370.)
      PARAMETER (PHINUL = 60. * RAD)     !MPA: Commented back in

c     MPA: Commented out below section
c      INCLUDE 'input.co'
c      IF(ibottom.eq.6)then
c         phinul = 30. * RAD
c      else
c         phinul = 60. * RAD
c      end if

C
C***FIRST EXECUTABLE STATEMENT SPH2GR.
C    
C  Error treatment.
C
      IFAIL = 0
      IF ((LATT .LE. -90.) .OR. (LATT .GT. 90.)) THEN
        IFAIL = 1
C       WRITE(*,*) '*** ERROR IN SPH2GR -- Illegal value of LATT ***'
        GOTO 999
      ENDIF
      IF ((DX .LE. 0.) .OR. (DY .LE. 0.)) THEN
        IFAIL = 2
        WRITE(*,*) '*** ERROR IN SPH2GR -- DX or DY not positive ***'
        write(*,*)DX,DY,LONG,LATT
        GOTO 999
      ENDIF
C
C  Convert degrees to radians.
C
      LAMBDA = LONG  * RAD
      PHI    = LATT  * RAD
      ALPHA  = YLONG * RAD
C
C  The computations.
C
      R = REARTH * COS(PHI) * (1 + SIN(PHINUL)) / (1 + SIN(PHI))
C
      XGRID = XPOLE + R * SIN(LAMBDA-ALPHA) / DX
      YGRID = YPOLE - R * COS(LAMBDA-ALPHA) / DY
C
  999 RETURN
C
C***END SPH2GR
C
      END subroutine sph2gr
