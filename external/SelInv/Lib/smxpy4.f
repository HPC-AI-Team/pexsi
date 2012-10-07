C***********************************************************************
C***********************************************************************
C
C   Version:        0.3
C   Last modified:  December 27, 1994
C   Authors:        Esmond G. Ng and Barry W. Peyton
C
C   Mathematical Sciences Section, Oak Ridge National Laboratory
C
C***********************************************************************
C***********************************************************************
C******     SMXPY4 .... MATRIX-VECTOR MULTIPLY            **************
C***********************************************************************
C***********************************************************************
C
C     PURPOSE - THIS ROUTINE PERFORMS A MATRIX-VECTOR MULTIPLY,
C               Y = Y + AX, ASSUMING DATA STRUCTURES USED IN
C               RECENTLY DEVELOPED SPARSE CHOLESKY CODES.  THE 
C               '4' SIGNIFIES LEVEL 4 LOOP UNROLLING.
C
C     INPUT PARAMETERS -
C        M      - NUMBER OF ROWS.
C        N      - NUMBER OF COLUMNS.
C        Y      - M-VECTOR TO WHICH AX WILL BE ADDED.
C        APNT   - INDEX VECTOR FOR A.  XA(I) POINTS TO THE
C                 FIRST NONZERO IN COLUMN I OF A.
C        Y      - ON OUTPUT, CONTAINS Y = Y + AX.
C
C***********************************************************************
C
      SUBROUTINE  SMXPY4 ( M, N, Y, APNT, A )
C
C***********************************************************************
C
C     -----------
C     PARAMETERS.
C     -----------
C
      INTEGER             M, N, LEVEL
C
      INTEGER             APNT(*)
C
      COMPLEX*16    Y(*), A(*)
C
      PARAMETER           ( LEVEL = 4 )
C
C     ----------------
C     LOCAL VARIABLES.
C     ----------------
C
      INTEGER             I, I1, I2, I3, I4,
     &                       J1, J2, J3, J4,
     &                    J, REMAIN
C
      COMPLEX*16    A1, A2, A3, A4
C
C***********************************************************************
C
      REMAIN = MOD ( N, LEVEL )
C
      GO TO ( 2000, 100, 200, 300 ), REMAIN+1
C
  100 CONTINUE
      J1 = APNT(1)
      I1 = APNT(1+1) - M
      A1 = - A(J1)*A(I1)
      DO  150  I = 1, M
          Y(I) = Y(I) + A1*A(I1)
          I1 = I1 + 1
  150 CONTINUE
      GO TO 2000
C
  200 CONTINUE
      J1 = APNT(1)
      J2 = APNT(1+1)
      I1 = J2        - M
      I2 = APNT(1+2) - M
      A1 = - A(J1)*A(I1)
      A2 = - A(J2)*A(I2)
      DO  250  I = 1, M
          Y(I) = ( (Y(I))
     &           + A1*A(I1)) + A2*A(I2)
          I1 = I1 + 1
          I2 = I2 + 1
  250 CONTINUE
      GO TO 2000
C
  300 CONTINUE
      J1 = APNT(1)
      J2 = APNT(1+1)
      J3 = APNT(1+2)
      I1 = J2        - M
      I2 = J3        - M
      I3 = APNT(1+3) - M
      A1 = - A(J1)*A(I1)
      A2 = - A(J2)*A(I2)
      A3 = - A(J3)*A(I3)
      DO  350  I = 1, M
          Y(I) = (( (Y(I))
     &           + A1*A(I1)) + A2*A(I2))
     &           + A3*A(I3)
          I1 = I1 + 1
          I2 = I2 + 1
          I3 = I3 + 1
  350 CONTINUE
      GO TO 2000
C
 2000 CONTINUE
      DO  4000  J = REMAIN+1, N, LEVEL
          J1 = APNT(J)
          J2 = APNT(J+1)
          J3 = APNT(J+2)
          J4 = APNT(J+3)
          I1 = J2        - M
          I2 = J3        - M
          I3 = J4        - M
          I4 = APNT(J+4) - M
          A1 = - A(J1)*A(I1)
          A2 = - A(J2)*A(I2)
          A3 = - A(J3)*A(I3)
          A4 = - A(J4)*A(I4)
          DO  3000  I = 1, M
              Y(I) = ((( (Y(I))
     &               + A1*A(I1)) + A2*A(I2))
     &               + A3*A(I3)) + A4*A(I4)
              I1 = I1 + 1
              I2 = I2 + 1
              I3 = I3 + 1
              I4 = I4 + 1
 3000     CONTINUE
 4000 CONTINUE
C
      RETURN
      END
