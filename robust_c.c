/*
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
*/
#define dmin(X,Y) ((X)<(Y)?(X):(Y))
#define dmax(X,Y) ((X)>(Y)?(X):(Y))
double hlqest_loc(const double *x, int *lb, int *rb, int *q, const int n)
{
    /* System generated locals */
    int i__1;
    double ret_val, r__1, r__2;

    /* Local variables */
    static int i__, j, k, l, k1, k2;
    static double am;
    static int iq, nn, sm, sq, lbi;
    static double amn;
    static int rbi;
    static double amx;
    extern double drand48(void);
    static int mdll, mdlu, ipiq;
    static int mdlrow;


/*     f2py intent(in) x */
/*     f2py intent(in) n                                                         C */

/*    REAL FUNCTION HLQEST */

/*    PURPOSE       COMPUTES THE HODGES-LEHMANN LOCATION ESTIMATOR: */
/*                  MEDIAN OF ( X(I) + X(J) ) / 2   FOR 1 LE I LE J LE N */

/*    USAGE         RESULT = HLQEST(X,N,LB,RB,Q) */

/*    ARGUMENTS  X   REAL ARRAY OF OBSERVATIONS  (INPUT) */

/*               N   INTEGER NUMBER OF OBSERVATIONS  (INPUT) */
/*                 * N MUST NOT BE LESS THAN 1 * */

/*   EXTERNAL ROUTINE */
/*              RANG  FUNCTION PROVIDING UNIFORM RANDOM VARIABLES */
/*                   IN THE INTERVAL (0,1) */
/*                   RANG REQUIRES A DUMMY INTEGER ARGUMENT */
/*              HSORT SORTS THE INPUT DATA */


/*   NOTES           HLQEST HAS AN EXPECTED TIME COMPLEXITY ON */
/*                   THE ORDER OF N * LG( N ) */

/*  J F MONAHAN, APRIL 1982, DEPT OF STAT, N C S U, RALEIGH, N C 27650 */
/*  FINAL VERSION  JUNE 1983 */



/*  TAKE CARE OF SPECIAL CASES: N=1 AND N=2 */

    /* Parameter adjustments */
    --q;
    --rb;
    --lb;
    --x;

    /* Function Body */
    if (n > 2) {
	goto L10;
    }
    ret_val = x[1];
    if (n == 1) {
	return ret_val;
    }
    ret_val = (x[1] + x[2]) / 2.f;
    return ret_val;

/*  FIND THE TOTAL NUMBER OF PAIRS (NN) AND THE MEDIAN(S) (K1,K2) NEEDED */

L10:
    nn = n * (n + 1) / 2;
    k1 = (nn + 1) / 2;
    k2 = (nn + 2) / 2;


/*  INITIALIZE LEFT AND RIGHT BOUNDS */

    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	lb[i__] = i__;
	rb[i__] = n;
/* L20: */
    }
/*  SM = NUMBER IN SET S AT STEP M */
    sm = nn;
/*  L = NUMBER OF PAIRS LESS THAN THOSE IN SET S AT STEP M */
    l = 0;


/*  USE THE MEDIAN OF X(I)'S TO PARTITION ON THE FIRST STEP */

    mdll = (n + 1) / 2;
    mdlu = (n + 2) / 2;
    am = x[mdll] + x[mdlu];
    goto L80;

/*  USE THE MIDRANGE OF SET S AS PARTITION ELEMENT WHEN TIES ARE LIKELY */
/*   -- OR GET THE AVERAGE OF THE LAST 2 ELEMENTS */

L30:
    amx = x[1] + x[1];
    amn = x[n] + x[n];
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*   SKIP THIS ROW IF NO ELEMENT IN IT IS IN SET S ON THIS STEP */
	if (lb[i__] > rb[i__]) {
	    goto L40;
	}
	lbi = lb[i__];
/*                             GET THE SMALLEST IN THIS ROW */
/* Computing MIN */
	r__1 = amn, r__2 = x[lbi] + x[i__];
	amn = dmin(r__1,r__2);
	rbi = rb[i__];
/*                             GET THE LARGEST IN THIS ROW */
/* Computing MAX */
	r__1 = amx, r__2 = x[rbi] + x[i__];
	amx = dmax(r__1,r__2);
L40:
	;
    }
    am = (amx + amn) / 2.f;
/*  BE CAREFUL TO CUT OFF SOMETHING -- ROUNDOFF CAN DO WIERD THINGS */
    if (am <= amn || am > amx) {
	am = amx;
    }
/*  UNLESS FINISHED, JUMP TO PARTITION STEP */
    if (amn != amx && sm != 2) {
	goto L80;
    }
/*  ALL DONE IF ALL OF S IS THE SAME -OR- IF ONLY 2 ELEMENTS ARE LEFT */
    ret_val = am / 2.f;
    return ret_val;

/*   *****   RESTART HERE UNLESS WORRIED ABOUT TIES   ***** */

L50:
/*                        USE RANDOM ROW MEDIAN AS PARTITION ELEMENT */
    k = (int) ((double) sm * drand48());
/*                        K IS A RANDOM INTEGER FROM O TO SM-1 */
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = i__;
	if (k <= rb[i__] - lb[i__]) {
	    goto L70;
	}
	k = k - rb[i__] + lb[i__] - 1;
/* L60: */
    }
/*                        J IS A RANDOM ROW --- NOW GET ITS MEDIAN */
L70:
    mdlrow = (lb[j] + rb[j]) / 2;
    am = x[j] + x[mdlrow];

/*       *****   PARTITION STEP   ***** */

/*  USE AM TO PARTITION S0 INTO 2 GROUPS: THOSE .LT. AM, THOSE .GE. AM */
/*  Q(I)= HOW MANY PAIRS (X(I)+X(J)) IN ROW I LESS THAN AM */
L80:
    j = n;
/*                              START IN UPPER RIGHT CORNER */
    sq = 0;
/*                              I COUNTS ROWS */
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	q[i__] = 0;
/*                              HAVE WE HIT THE DIAGONAL ? */
L90:
	if (j < i__) {
	    goto L110;
	}
/*                              SHALL WE MOVE LEFT ? */
	if (x[i__] + x[j] < am) {
	    goto L100;
	}
	--j;
	goto L90;
/*                              WE'RE DONE IN THIS ROW */
L100:
	q[i__] = j - i__ + 1;
/*  SQ = TOTAL NUMBER OF PAIRS LESS THAN AM */
	sq += q[i__];
L110:
	;
    }

/*  ***  FINISHED PARTITION --- START BRANCHING  *** */

/*  IF CONSECUTIVE PARTITIONS ARE THE SAME WE PROBABLY HAVE TIES */
    if (sq == l) {
	goto L30;
    }

/*  ARE WE NEARLY DONE, WITH THE VALUES WE WANT ON THE BORDER? */
/*  IF(WE NEED  MAX OF THOSE .LT. AM -OR- MIN OF THOSE .GE. AM) GO TO 90 */

    if (sq == k2 - 1) {
	goto L180;
    }

/*  THE SET S IS SPLIT, WHICH PIECE DO WE KEEP? */
/*  70  =  CUT OFF BOTTOM,   90  =  NEARLY DONE,   60  =  CUT OFF TOP */

    if ((i__1 = sq - k1) < 0) {
	goto L140;
    } else if (i__1 == 0) {
	goto L180;
    } else {
	goto L120;
    }

/*  NEW S = (OLD S) .INTERSECT. (THOSE .LT. AM) */
L120:
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*                            RESET RIGHT BOUNDS FOR EACH ROW */
	rb[i__] = i__ + q[i__] - 1;
/* L130: */
    }
    goto L160;
/*  NEW S = (OLD S) .INTERSECT. (THOSE .GE. AM) */
L140:
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*                            RESET LEFT BOUNDS FOR EACH ROW */
	lb[i__] = i__ + q[i__];
/* L150: */
    }

/*  COUNT   SM = NUMBER OF PAIRS STILL IN NEW SET S */
/*           L = NUMBER OF PAIRS LESS THAN THOSE IN NEW SET S */
L160:
    l = 0;
    sm = 0;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l = l + lb[i__] - i__;
	sm = sm + rb[i__] - lb[i__] + 1;
/* L170: */
    }

/*        *****   NORMAL RESTART JUMP   ***** */

    if (sm > 2) {
	goto L50;
    }
/*  CAN ONLY GET TO 2 LEFT IF K1.NE.K2  -- GO GET THEIR AVERAGE */
    goto L30;

/*  FIND:   MAX OF THOSE .LT. AM */
/*          MIN OF THOSE .GE. AM */
L180:
    amn = x[n] + x[n];
    amx = x[1] + x[1];
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iq = q[i__];
	ipiq = i__ + iq;
	if (iq > 0) {
/* Computing MAX */
	    r__1 = amx, r__2 = x[i__] + x[ipiq - 1];
	    amx = dmax(r__1,r__2);
	}
	ipiq = i__ + iq;
	if (iq < n - i__ + 1) {
/* Computing MIN */
	    r__1 = amn, r__2 = x[i__] + x[ipiq];
	    amn = dmin(r__1,r__2);
	}
/* L190: */
    }
    ret_val = (amn + amx) / 4.f;
/*  WE ARE DONE, BUT WHICH SITUATION ARE WE IN? */
    if (k1 < k2) {
	return ret_val;
    }
    if (sq == k1) {
	ret_val = amx / 2.f;
    }
    if (sq == k1 - 1) {
	ret_val = amn / 2.f;
    }
    return ret_val;
} /* hlqest_ loc */

double hlqest(const double* x, const size_t n) {
  double retval, *y;
  int *lb, *rb, *q;
  if ( (lb=(int *) malloc(sizeof(int)*n)) == NULL) {
    printf("Failed to allocate lb in hlqest()\n");
    return 0.0/0.0;
  }
  if ( (rb=(int *) malloc(sizeof(int)*n)) == NULL) {
    printf("Failed to allocate rb in hlqest()\n");
    return 0.0/0.0;
  }
  if ( (q=(int *) malloc(sizeof(int)*n)) == NULL) {
    printf("Failed to allocate q in hlqest()\n");
    return 0.0/0.0;
  }
  if ( (y=(double *) malloc(sizeof(double)*n)) == NULL) {
    printf("Failed to allocate y in hlqest()\n");
    return 0.0/0.0;
  }


  /* back up the list y */
  memcpy(y, x, n*sizeof(double));

  /* hlqest_loc requires a sorted list */
  qsort(y, n, sizeof(double), cmp_double);
  retval=hlqest_loc(y, lb, rb, q, n);

  free((void *) lb);
  free((void *) rb);
  free((void *) q);
  return(retval);
 
}
