      SUBROUTINE VECG2S(XGRID, YGRID, UGRID, VGRID, XPOLE, YPOLE,
     >                  USPH, VSPH, IFAIL)
C
C***BEGIN PROLOGUE VECG2S
C***DATE WRITTEN    4 FEB. 1991
C***REVISION DATE   4 FEB. 1991
C***AUTHOR
C            Bj{\o}rn {\AA}dlandsvik
C            Institute of Marine Research,
C            Postboks 1870 Nordnes
C            N-5024 Bergen, Norway
C            Email..  bjorn@imr.no
C
C***PURPOSE
C  Computes the east and north components of a current vector
C  given by x-y components in the directions of a polar stereographic 
C  grid coordinate system.
C
C***PARAMETERS
C***ON ENTRY
C     XGRID  Real,
C        Grid coordinate in X-direction of the vectors footpoint.
C     YGRID  Real,
C        Grid coordinate in Y-direction of the vectors footpoint.
C     UGRID  Real,
C        Grid component of vector -- X direction.
C        Arbitrary units are allowed.
C     VGRID  Real,
C        Grid component of vector -- Y direction.
C        In same unit as UGRID.
C     XPOLE  Real,	
C  	 Grid coordinate of North Pole in X-direction.
C     YPOLE  Real,	
C  	 Grid coordinate of North Pole in Y-direction.
C
C***ON RETURN
C     USPH   Real,
C        Vector component -- Easterly direction.
C        Given in same unit as UGRID and VGRID.
C     VSPH   Real,
C        Vector component -- Northerly direction.    
C        Given in same unit as UGRID and VGRID.
C     IFAIL  Integer,
C        Error indicator. 
C        Returns : 1 if (XGRID, YGRID) = North Pole,
C                  0 otherwise.
C
C***DESCRIPTION
C  The grid coordinate system is locally rotated an angle a,
C  thus the velocity vector is transformed by the matrix formula.
C
C    | USPH |    | cos(a)  -sin(a) |   | UGRID |
C    | VSPH | =  | sin(a)   cos(a) | * | VGRID | .
C  
C  The angle a is given by:
C     tan(a) = (XPOLE-XGRID) / (YPOLE-YGRID) .
C 
C  Using R**2 = (XPOLE-XGRID)**2 + (YPOLE-YGRID)**2
C  the following expressions are used:
C
C  sin(a) = (XPOLE-XGRID)/R, cos(a) = (YPOLE-YGRID)/R . 
C
C  Restriction: 
C    The vector's footpoint can not be the North Pole,
C    because East and North directions are not defined there.
C 
C***LOCAL VARIABLES
C     R		Real,
C        Distance from vectors footpoint to North Pole.
C     SINA      Real,
C        Sine of rotation angle.
C     COSA	Real,
C        Cosine of rotation angle.
C
C***RELATED ROUTINES
C  The vector's footpoint can be transformed with GR2SPH.
C  This routine is the inverse of VECS2G.
C
C***ROUTINES CALLED-NONE
C***END PROLOGUE VECG2S
C
C  Parameters
C
      REAL XGRID, YGRID, UGRID, VGRID
      REAL XPOLE, YPOLE
      REAL USPH, VSPH
      INTEGER IFAIL
C
C  Local variables.
C
      REAL R, SINA, COSA
C
C***FIRST EXECUTABLE STATEMENT VECG2S
C    
C  Error treatment.
C
      IFAIL = 0
      IF ((XGRID .EQ. XPOLE) .AND. (YGRID .EQ. YPOLE)) THEN
        IFAIL = 1
        WRITE(*,*) '*** ERROR IN VECG2S -- North Pole is illegal ***'
        GOTO 999
      ENDIF
C
C  Compute R, SINA, COSA
C
      R = SQRT((XPOLE-XGRID)*(XPOLE-XGRID)+(YPOLE-YGRID)*(YPOLE-YGRID))
      SINA = (XPOLE-XGRID) / R
      COSA = (YPOLE-YGRID) / R
C
C  Coordinate transformation.
C
      USPH = COSA * UGRID - SINA * VGRID
      VSPH = SINA * UGRID + COSA * VGRID
C
 999  RETURN
C
C***END VECG2S
C
      END
