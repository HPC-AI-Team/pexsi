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
C******     SMXPY8 .... MATRIX-VECTOR MULTIPLY            **************
C***********************************************************************
C***********************************************************************
C
C     PURPOSE - THIS ROUTINE PERFORMS A MATRIX-VECTOR MULTIPLY,
C               Y = Y + AX, ASSUMING DATA STRUCTURES USED IN
C               RECENTLY DEVELOPED SPARSE CHOLESKY CODES.  THE 
C               '8' SIGNIFIES LEVEL 8 LOOP UNROLLING.
C
C     INPUT PARAMETERS -
C        M      - NUMBER OF ROWS.
C        N      - NUMBER OF COLUMNS.
C        Y      - M-VECTOR TO WHICH AX WILL BE ADDED.
C        APNT   - INDEX VECTOR FOR A.  APNT(I) POINTS TO THE
C                 FIRST NONZERO IN COLUMN I OF A.
C        Y      - ON OUTPUT, CONTAINS Y = Y + AX.
C
C***********************************************************************
C
      SUBROUTINE  SMXPY8 ( M, N, Y, APNT, A )
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
      PARAMETER           ( LEVEL = 8 )
C
C     ----------------
C     LOCAL VARIABLES.
C     ----------------
C
      INTEGER             I, I1, I2, I3, I4, I5, I6, I7, I8,
     &                       J1, J2, J3, J4, J5, J6, J7, J8,
     &                    J, REMAIN
C
      COMPLEX*16    A1, A2, A3, A4, A5, A6, A7, A8
C
C***********************************************************************
C
      REMAIN = MOD ( N, LEVEL )
C
      GO TO ( 2000, 100, 200, 300,
     &         400, 500, 600, 700  ), REMAIN+1
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
  400 CONTINUE
      J1 = APNT(1)
      J2 = APNT(1+1)
      J3 = APNT(1+2)
      J4 = APNT(1+3)
      I1 = J2        - M
      I2 = J3        - M
      I3 = J4        - M
      I4 = APNT(1+4) - M
      A1 = - A(J1)*A(I1)
      A2 = - A(J2)*A(I2)
      A3 = - A(J3)*A(I3)
      A4 = - A(J4)*A(I4)
      DO  450  I = 1, M
          Y(I) = ((( (Y(I))
     &           + A1*A(I1)) + A2*A(I2))
     &           + A3*A(I3)) + A4*A(I4)
          I1 = I1 + 1
          I2 = I2 + 1
          I3 = I3 + 1
          I4 = I4 + 1
  450 CONTINUE
      GO TO 2000
C
  500 CONTINUE
      J1 = APNT(1)
      J2 = APNT(1+1)
      J3 = APNT(1+2)
      J4 = APNT(1+3)
      J5 = APNT(1+4)
      I1 = J2        - M
      I2 = J3        - M
      I3 = J4        - M
      I4 = J5        - M
      I5 = APNT(1+5) - M
      A1 = - A(J1)*A(I1)
      A2 = - A(J2)*A(I2)
      A3 = - A(J3)*A(I3)
      A4 = - A(J4)*A(I4)
      A5 = - A(J5)*A(I5)
      DO  550  I = 1, M
          Y(I) = (((( (Y(I))
     &           + A1*A(I1)) + A2*A(I2))
     &           + A3*A(I3)) + A4*A(I4))
     &           + A5*A(I5)
          I1 = I1 + 1
          I2 = I2 + 1
          I3 = I3 + 1
          I4 = I4 + 1
          I5 = I5 + 1
  550 CONTINUE
      GO TO 2000
C
  600 CONTINUE
      J1 = APNT(1)
      J2 = APNT(1+1)
      J3 = APNT(1+2)
      J4 = APNT(1+3)
      J5 = APNT(1+4)
      J6 = APNT(1+5)
      I1 = J2        - M
      I2 = J3        - M
      I3 = J4        - M
      I4 = J5        - M
      I5 = J6        - M
      I6 = APNT(1+6) - M
      A1 = - A(J1)*A(I1)
      A2 = - A(J2)*A(I2)
      A3 = - A(J3)*A(I3)
      A4 = - A(J4)*A(I4)
      A5 = - A(J5)*A(I5)
      A6 = - A(J6)*A(I6)
      DO  650  I = 1, M
          Y(I) = ((((( (Y(I))
     &           + A1*A(I1)) + A2*A(I2))
     &           + A3*A(I3)) + A4*A(I4))
     &           + A5*A(I5)) + A6*A(I6)
          I1 = I1 + 1
          I2 = I2 + 1
          I3 = I3 + 1
          I4 = I4 + 1
          I5 = I5 + 1
          I6 = I6 + 1
  650 CONTINUE
      GO TO 2000
C
  700 CONTINUE
      J1 = APNT(1)
      J2 = APNT(1+1)
      J3 = APNT(1+2)
      J4 = APNT(1+3)
      J5 = APNT(1+4)
      J6 = APNT(1+5)
      J7 = APNT(1+6)
      I1 = J2        - M
      I2 = J3        - M
      I3 = J4        - M
      I4 = J5        - M
      I5 = J6        - M
      I6 = J7        - M
      I7 = APNT(1+7) - M
      A1 = - A(J1)*A(I1)
      A2 = - A(J2)*A(I2)
      A3 = - A(J3)*A(I3)
      A4 = - A(J4)*A(I4)
      A5 = - A(J5)*A(I5)
      A6 = - A(J6)*A(I6)
      A7 = - A(J7)*A(I7)
      DO  750  I = 1, M
          Y(I) = (((((( (Y(I))
     &           + A1*A(I1)) + A2*A(I2))
     &           + A3*A(I3)) + A4*A(I4))
     &           + A5*A(I5)) + A6*A(I6))
     &           + A7*A(I7)
          I1 = I1 + 1
          I2 = I2 + 1
          I3 = I3 + 1
          I4 = I4 + 1
          I5 = I5 + 1
          I6 = I6 + 1
          I7 = I7 + 1
  750 CONTINUE
      GO TO 2000
C
 2000 CONTINUE
      DO  4000  J = REMAIN+1, N, LEVEL
          J1 = APNT(J)
          J2 = APNT(J+1)
          J3 = APNT(J+2)
          J4 = APNT(J+3)
          J5 = APNT(J+4)
          J6 = APNT(J+5)
          J7 = APNT(J+6)
          J8 = APNT(J+7)
          I1 = J2        - M
          I2 = J3        - M
          I3 = J4        - M
          I4 = J5        - M
          I5 = J6        - M
          I6 = J7        - M
          I7 = J8        - M
          I8 = APNT(J+8) - M
          A1 = - A(J1)*A(I1)
          A2 = - A(J2)*A(I2)
          A3 = - A(J3)*A(I3)
          A4 = - A(J4)*A(I4)
          A5 = - A(J5)*A(I5)
          A6 = - A(J6)*A(I6)
          A7 = - A(J7)*A(I7)
          A8 = - A(J8)*A(I8)
          DO  3000  I = 1, M
              Y(I) = ((((((( (Y(I))
     &               + A1*A(I1)) + A2*A(I2))
     &               + A3*A(I3)) + A4*A(I4))
     &               + A5*A(I5)) + A6*A(I6))
     &               + A7*A(I7)) + A8*A(I8)
              I1 = I1 + 1
              I2 = I2 + 1
              I3 = I3 + 1
              I4 = I4 + 1
              I5 = I5 + 1
              I6 = I6 + 1
              I7 = I7 + 1
              I8 = I8 + 1
 3000     CONTINUE
 4000 CONTINUE
C
      RETURN
      END
