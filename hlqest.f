 C     FILE: hlqest.f
C     ALGORITHM 616 COLLECTED ALGORITHMS FROM ACM.
C     ALGORITHM APPEARED IN ACM-TRANS. MATH. SOFTWARE, VOL.10, NO. 3,
C     SEP., 1984, P. 265-270.
C
C     HLQEST: a Fortran subprogram for computing the Hodges-Lehmann
C     location estimator, median of ( x(i) + x(j) )/2 for
C     1 .LE. I .LE. J .LE. n.
C     (See J.F. Monahan, ACM TOMS 10 (1984) 265-270).
C
C      
      SUBROUTINE HSORT(K, N)                                           
C     f2py intent(inout) k                                                      C     f2py intent(in) n                                                              
C  HEAPSORT ALGORITHM FOR SORTING ON VECTOR OF KEYS K OF LENGTH N
C  J F MONAHAN        TRANSCRIBED FROM KNUTH, VOL 2, PP 146-7.
      REAL K(N), KK
      INTEGER R
      IF (N.LE.1) RETURN
      L = N/2 + 1
      R = N
 10    IF (L.GT.1) GO TO 20
      KK = K(R)
      K(R) = K(1)
      R = R - 1
      IF (R.EQ.1) GO TO 80
      GO TO 30
 20    L = L - 1
      KK = K(L)
 30    J = L
 40     I = J
      J = 2*J
      IF (J-R) 50, 60, 70
 50    IF (K(J).LT.K(J+1)) J = J + 1
 60     IF (KK.GT.K(J)) GO TO 70
      K(I) = K(J)
      GO TO 40
 70    K(I) = KK
      GO TO 10
 80    K(1) = KK
      RETURN
      END
C
      REAL FUNCTION RANG(IXX)                                          
C  UNIFORM PSEUDORANDOM NUMBER GENERATOR
C  FORTRAN VERSION OF LEWIS, GOODMAN, MILLER
C  SCHRAGE,  ACM TOMS V.5 (1979) P132
C  FIRST CALL SETS SEED TO IXX, LATER IXX IGNORED
C     f2py intent(inout) ixx                                                              
      INTEGER A, P, IX, B15, B16, XHI, XALO, LEFTLO, FHI, K
      DATA A /16807/, B15 /32768/, B16 /65536/, P /2147483647/
      DATA IX /0/
      IF (IX.EQ.0) IX = IXX
      XHI = IX/B16
      XALO = (IX-XHI*B16)*A
      LEFTLO = XALO/B16
      FHI = XHI*A + LEFTLO
      K = FHI/B15
      IX = (((XALO-LEFTLO*B16)-P)+(FHI-K*B15)*B16) + K
      IF (IX.LT.0) IX = IX + P
      RANG = FLOAT(IX)*4.656612875E-10
      RETURN
      END

      REAL FUNCTION HLQEST(X, N)
C     f2py intent(in) x
C     f2py intent(in) n                                                         C
C
C    REAL FUNCTION HLQEST
C
C    PURPOSE       COMPUTES THE HODGES-LEHMANN LOCATION ESTIMATOR:
C                  MEDIAN OF ( X(I) + X(J) ) / 2   FOR 1 LE I LE J LE N
C
C    USAGE         RESULT = HLQEST(X,N,LB,RB,Q)
C
C    ARGUMENTS  X   REAL ARRAY OF OBSERVATIONS  (INPUT)
C
C               N   INTEGER NUMBER OF OBSERVATIONS  (INPUT)
C                 * N MUST NOT BE LESS THAN 1 *
C
C   EXTERNAL ROUTINE
C              RANG  FUNCTION PROVIDING UNIFORM RANDOM VARIABLES
C                   IN THE INTERVAL (0,1)
C                   RANG REQUIRES A DUMMY INTEGER ARGUMENT
C              HSORT SORTS THE INPUT DATA
C
C
C   NOTES           HLQEST HAS AN EXPECTED TIME COMPLEXITY ON
C                   THE ORDER OF N * LG( N )
C
C  J F MONAHAN, APRIL 1982, DEPT OF STAT, N C S U, RALEIGH, N C 27650
C  FINAL VERSION  JUNE 1983
C
C      
      REAL X(N), AMN, AMX, AM
      INTEGER LB(N), RB(N), Q(N), SM, SQ, I, J, K, K1, K2, L, N, NN,
     * MDLL, MDLU, LBI, RBI, MDLROW, IQ
C
C  TAKE CARE OF SPECIAL CASES: N=1 AND N=2
C
      IF (N.GT.2) GO TO 10
      HLQEST = X(1)
      IF (N.EQ.1) RETURN
      HLQEST = (X(1)+X(2))/2.
      RETURN
C
C  FIND THE TOTAL NUMBER OF PAIRS (NN) AND THE MEDIAN(S) (K1,K2) NEEDED
C
   10 NN = (N*(N+1))/2
      K1 = (NN+1)/2
      K2 = (NN+2)/2
C
C  SORT THE DATA
C
      CALL HSORT(X,N)
C
C  INITIALIZE LEFT AND RIGHT BOUNDS
C
      DO 20 I=1,N
        LB(I) = I
        RB(I) = N
   20 CONTINUE
C  SM = NUMBER IN SET S AT STEP M
      SM = NN
C  L = NUMBER OF PAIRS LESS THAN THOSE IN SET S AT STEP M
      L = 0
C
C
C  USE THE MEDIAN OF X(I)'S TO PARTITION ON THE FIRST STEP
C
      MDLL = (N+1)/2
      MDLU = (N+2)/2
      AM = X(MDLL) + X(MDLU)
      GO TO 80
C
C  USE THE MIDRANGE OF SET S AS PARTITION ELEMENT WHEN TIES ARE LIKELY
C   -- OR GET THE AVERAGE OF THE LAST 2 ELEMENTS
C
   30 AMX = X(1) + X(1)
      AMN = X(N) + X(N)
      DO 40 I=1,N
C   SKIP THIS ROW IF NO ELEMENT IN IT IS IN SET S ON THIS STEP
        IF (LB(I).GT.RB(I)) GO TO 40
        LBI = LB(I)
C                             GET THE SMALLEST IN THIS ROW
        AMN = AMIN1(AMN,X(LBI)+X(I))
        RBI = RB(I)
C                             GET THE LARGEST IN THIS ROW
        AMX = AMAX1(AMX,X(RBI)+X(I))
   40 CONTINUE
      AM = (AMX+AMN)/2.
C  BE CAREFUL TO CUT OFF SOMETHING -- ROUNDOFF CAN DO WIERD THINGS
      IF (AM.LE.AMN .OR. AM.GT.AMX) AM = AMX
C  UNLESS FINISHED, JUMP TO PARTITION STEP
      IF (AMN.NE.AMX .AND. SM.NE.2) GO TO 80
C  ALL DONE IF ALL OF S IS THE SAME -OR- IF ONLY 2 ELEMENTS ARE LEFT
      HLQEST = AM/2.
      RETURN
C
C   *****   RESTART HERE UNLESS WORRIED ABOUT TIES   *****
C
   50 CONTINUE
C                        USE RANDOM ROW MEDIAN AS PARTITION ELEMENT
      K = IFIX(FLOAT(SM)*RANG(SM))
C                        K IS A RANDOM INTEGER FROM O TO SM-1
      DO 60 I=1,N
        J = I
        IF (K.LE.RB(I)-LB(I)) GO TO 70
        K = K - RB(I) + LB(I) - 1
   60 CONTINUE
C                        J IS A RANDOM ROW --- NOW GET ITS MEDIAN
   70 MDLROW = (LB(J)+RB(J))/2
      AM = X(J) + X(MDLROW)
C
C       *****   PARTITION STEP   *****
C
C  USE AM TO PARTITION S0 INTO 2 GROUPS: THOSE .LT. AM, THOSE .GE. AM
C  Q(I)= HOW MANY PAIRS (X(I)+X(J)) IN ROW I LESS THAN AM
   80 CONTINUE
      J = N
C                              START IN UPPER RIGHT CORNER
      SQ = 0
C                              I COUNTS ROWS
      DO 110 I=1,N
        Q(I) = 0
C                              HAVE WE HIT THE DIAGONAL ?
   90   IF (J.LT.I) GO TO 110
C                              SHALL WE MOVE LEFT ?
        IF (X(I)+X(J).LT.AM) GO TO 100
        J = J - 1
        GO TO 90
C                              WE'RE DONE IN THIS ROW
  100   Q(I) = J - I + 1
C  SQ = TOTAL NUMBER OF PAIRS LESS THAN AM
        SQ = SQ + Q(I)
  110 CONTINUE
C
C  ***  FINISHED PARTITION --- START BRANCHING  ***
C
C  IF CONSECUTIVE PARTITIONS ARE THE SAME WE PROBABLY HAVE TIES
      IF (SQ.EQ.L) GO TO 30
C
C  ARE WE NEARLY DONE, WITH THE VALUES WE WANT ON THE BORDER?
C  IF(WE NEED  MAX OF THOSE .LT. AM -OR- MIN OF THOSE .GE. AM) GO TO 90
C
      IF (SQ.EQ.K2-1) GO TO 180
C
C  THE SET S IS SPLIT, WHICH PIECE DO WE KEEP?
C  70  =  CUT OFF BOTTOM,   90  =  NEARLY DONE,   60  =  CUT OFF TOP
C
      IF (SQ-K1) 140, 180, 120
C
C  NEW S = (OLD S) .INTERSECT. (THOSE .LT. AM)
  120 CONTINUE
      DO 130 I=1,N
C                            RESET RIGHT BOUNDS FOR EACH ROW
        RB(I) = I + Q(I) - 1
  130 CONTINUE
      GO TO 160
C  NEW S = (OLD S) .INTERSECT. (THOSE .GE. AM)
  140 CONTINUE
      DO 150 I=1,N
C                            RESET LEFT BOUNDS FOR EACH ROW
        LB(I) = I + Q(I)
  150 CONTINUE
C
C  COUNT   SM = NUMBER OF PAIRS STILL IN NEW SET S
C           L = NUMBER OF PAIRS LESS THAN THOSE IN NEW SET S
  160 L = 0
      SM = 0
      DO 170 I=1,N
        L = L + LB(I) - I
        SM = SM + RB(I) - LB(I) + 1
  170 CONTINUE
C
C        *****   NORMAL RESTART JUMP   *****
C
      IF (SM.GT.2) GO TO 50
C  CAN ONLY GET TO 2 LEFT IF K1.NE.K2  -- GO GET THEIR AVERAGE
      GO TO 30
C
C  FIND:   MAX OF THOSE .LT. AM
C          MIN OF THOSE .GE. AM
  180 CONTINUE
      AMN = X(N) + X(N)
      AMX = X(1) + X(1)
      DO 190 I=1,N
        IQ = Q(I)
        IPIQ = I + IQ
        IF (IQ.GT.0) AMX = AMAX1(AMX,X(I)+X(IPIQ-1))
        IPIQ = I + IQ
        IF (IQ.LT.N-I+1) AMN = AMIN1(AMN,X(I)+X(IPIQ))
  190 CONTINUE
      HLQEST = (AMN+AMX)/4.
C  WE ARE DONE, BUT WHICH SITUATION ARE WE IN?
      IF (K1.LT.K2) RETURN
      IF (SQ.EQ.K1) HLQEST = AMX/2.
      IF (SQ.EQ.K1-1) HLQEST = AMN/2.
      RETURN
      END

      PROGRAM TEST
      REAL XX(1000)
      INTEGER I
      ID=1000
      DO I=1,1000
         XX(I)=1000-I
      END DO
C     XX(1000)=500
      PRINT *,HLQEST(XX,1000)
      END
