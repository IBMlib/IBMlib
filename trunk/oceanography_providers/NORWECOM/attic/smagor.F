#include "cppdefs.h"
      SUBROUTINE SMAGOR(IM,JM,KB,AAM,AAM2D,HORCON,DX,DY,
     +                  U,V,DZ)
C***BEGIN PROLOGUE SMAGOR
C***DATE WRITTEN   911010 (YYMMDD)
C***REVISION DATE  920129 (YYMMDD)
C   In the last revision the diffusivities are extended
C   to boundary cells.
C***AUTHOR
C            Jarle Berntsen, Institute of Marine Research,
C            Postboks 1870,
C            N-5024 Bergen-Nordnes, Norway
C            Email..  jarleb@imr.no
C
C***PURPOSE SMAGOR computes the horizontal diffusivities
C           AAM and AAM2D with Smagorinsky's technique.
C
C***ROUTINES CALLED-NONE
C***END PROLOGUE SMAGOR
C
C   Global variables.
C
      INTEGER IM,JM,KB
      REAL AAM(IM,JM,KB),AAM2D(IM,JM),HORCON,DX,DY,
     +     U(IM,JM,KB),V(IM,JM,KB),DZ(KB)
C
C   Local variables.
C
      INTEGER I,J,K
C
      DO 10 K = 1,KB-1
        DO 10 J = 2,JM-1
          DO 10 I = 2,IM-1
            AAM(I,J,K)=HORCON*DX*DY
     1  *SQRT( ((U(I+1,J,K)-U(I,J,K))/DX)**2
     2        +((V(I,J+1,K)-V(I,J,K))/DY)**2
     3+.5*(.25*(U(I,J+1,K)+U(I+1,J+1,K)-U(I,J-1,K)-U(I+1,J-1,K))
     4     /DY   
     5     + .25*(V(I+1,J,K)+V(I+1,J+1,K)-V(I-1,J,K)-V(I-1,J+1,K))
     6     /DX)**2)
10    CONTINUE
C
C  Extend the diffusivities to the boundary grid cells.
C
      DO 20 K = 1,KB-1
        DO 20 J = 2,JM-1
          AAM(1,J,K) = AAM(2,J,K)
          AAM(IM,J,K) = AAM(IM-1,J,K)
20    CONTINUE
      DO 30 K = 1,KB-1
        DO 30 I = 2,IM-1
          AAM(I,1,K) = AAM(I,2,K)
          AAM(I,JM,K) = AAM(I,JM-1,K)
30    CONTINUE
      DO 40 K = 1,KB-1
        AAM(1,1,K) = AAM(2,2,K)
        AAM(1,JM,K) = AAM(2,JM-1,K)
        AAM(IM,JM,K) = AAM(IM-1,JM-1,K)
        AAM(IM,1,K) = AAM(IM-1,2,K)
40    CONTINUE
C
      DO 100 J=1,JM
        DO 100 I=1,IM
          AAM2D(I,J)=0.
100   CONTINUE
      DO 110 K = 1,KB-1
        DO 110 J = 1,JM
          DO 110 I = 1,IM
            AAM2D(I,J)=AAM2D(I,J)+AAM(I,J,K)*DZ(K)
110   CONTINUE
C
      RETURN
      END
