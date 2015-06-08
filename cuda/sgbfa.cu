#ifdef __cplusplus
extern "C" {
#endif

#include "f2c.h"

/* DECK SGBFA */
/* Subroutine */
__device__ int sgbfa_(real *abd, integer *lda, integer *n, integer *ml, 
	integer *mu, integer *ipvt, integer *info)
{
    /* Table of constant values */

    integer c__1 = 1;

    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    integer i__, j, k, l, m;
    real t;
    integer i0, j0, j1, lm, mm, ju, jz, kp1, nm1;
    extern /* Subroutine */ __device__ int sscal_(integer *, real *, real *, integer *), 
	    saxpy_(integer *, real *, real *, integer *, real *, integer *);
    extern __device__ integer isamax_(integer *, real *, integer *);

/* ***BEGIN PROLOGUE  SGBFA */
/* ***PURPOSE  Factor a band matrix using Gaussian elimination. */
/* ***CATEGORY  D2A2 */
/* ***TYPE      SINGLE PRECISION (SGBFA-S, DGBFA-D, CGBFA-C) */
/* ***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION */
/* ***AUTHOR  Moler, C. B., (U. of New Mexico) */
/* ***DESCRIPTION */

/*     SGBFA factors a real band matrix by elimination. */

/*     SGBFA is usually called by SBGCO, but it can be called */
/*     directly with a saving in time if  RCOND  is not needed. */

/*     On Entry */

/*        ABD     REAL(LDA, N) */
/*                contains the matrix in band storage.  The columns */
/*                of the matrix are stored in the columns of  ABD  and */
/*                the diagonals of the matrix are stored in rows */
/*                ML+1 through 2*ML+MU+1 of  ABD . */
/*                See the comments below for details. */

/*        LDA     INTEGER */
/*                the leading dimension of the array  ABD . */
/*                LDA must be .GE. 2*ML + MU + 1 . */

/*        N       INTEGER */
/*                the order of the original matrix. */

/*        ML      INTEGER */
/*                number of diagonals below the main diagonal. */
/*                0 .LE. ML .LT. N . */

/*        MU      INTEGER */
/*                number of diagonals above the main diagonal. */
/*                0 .LE. MU .LT. N . */
/*                More efficient if  ML .LE. MU . */
/*     On Return */

/*        ABD     an upper triangular matrix in band storage and */
/*                the multipliers which were used to obtain it. */
/*                The factorization can be written  A = L*U , where */
/*                L  is a product of permutation and unit lower */
/*                triangular matrices and  U  is upper triangular. */

/*        IPVT    INTEGER(N) */
/*                an integer vector of pivot indices. */

/*        INFO    INTEGER */
/*                = 0  normal value. */
/*                = K  if  U(K,K) .EQ. 0.0 .  This is not an error */
/*                     condition for this subroutine, but it does */
/*                     indicate that SGBSL will divide by zero if */
/*                     called.  Use  RCOND  in SBGCO for a reliable */
/*                     indication of singularity. */

/*     Band Storage */

/*           If  A  is a band matrix, the following program segment */
/*           will set up the input. */

/*                   ML = (band width below the diagonal) */
/*                   MU = (band width above the diagonal) */
/*                   M = ML + MU + 1 */
/*                   DO 20 J = 1, N */
/*                      I1 = MAX(1, J-MU) */
/*                      I2 = MIN(N, J+ML) */
/*                      DO 10 I = I1, I2 */
/*                         K = I - J + M */
/*                         ABD(K,J) = A(I,J) */
/*                10    CONTINUE */
/*                20 CONTINUE */

/*           This uses rows  ML+1  through  2*ML+MU+1  of  ABD . */
/*           In addition, the first  ML  rows in  ABD  are used for */
/*           elements generated during the triangularization. */
/*           The total number of rows needed in  ABD  is  2*ML+MU+1 . */
/*           The  ML+MU by ML+MU  upper left triangle and the */
/*           ML by ML  lower right triangle are not referenced. */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  ISAMAX, SAXPY, SSCAL */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  SGBFA */


/* ***FIRST EXECUTABLE STATEMENT  SGBFA */
    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    --ipvt;

    /* Function Body */
    m = *ml + *mu + 1;
    *info = 0;

/*     ZERO INITIAL FILL-IN COLUMNS */

    j0 = *mu + 2;
    j1 = min(*n,m) - 1;
    if (j1 < j0) {
	goto L30;
    }
    i__1 = j1;
    for (jz = j0; jz <= i__1; ++jz) {
	i0 = m + 1 - jz;
	i__2 = *ml;
	for (i__ = i0; i__ <= i__2; ++i__) {
	    abd[i__ + jz * abd_dim1] = 0.f;
/* L10: */
	}
/* L20: */
    }
L30:
    jz = j1;
    ju = 0;

/*     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING */

    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L130;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;

/*        ZERO NEXT FILL-IN COLUMN */

	++jz;
	if (jz > *n) {
	    goto L50;
	}
	if (*ml < 1) {
	    goto L50;
	}
	i__2 = *ml;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    abd[i__ + jz * abd_dim1] = 0.f;
/* L40: */
	}
L50:

/*        FIND L = PIVOT INDEX */

/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	lm = min(i__2,i__3);
	i__2 = lm + 1;
	l = isamax_(&i__2, &abd[m + k * abd_dim1], &c__1) + m - 1;
	ipvt[k] = l + k - m;

/*        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED */

	if (abd[l + k * abd_dim1] == 0.f) {
	    goto L100;
	}

/*           INTERCHANGE IF NECESSARY */

	if (l == m) {
	    goto L60;
	}
	t = abd[l + k * abd_dim1];
	abd[l + k * abd_dim1] = abd[m + k * abd_dim1];
	abd[m + k * abd_dim1] = t;
L60:

/*           COMPUTE MULTIPLIERS */

	t = -1.f / abd[m + k * abd_dim1];
	sscal_(&lm, &t, &abd[m + 1 + k * abd_dim1], &c__1);

/*           ROW ELIMINATION WITH COLUMN INDEXING */

/* Computing MIN */
/* Computing MAX */
	i__3 = ju, i__4 = *mu + ipvt[k];
	i__2 = max(i__3,i__4);
	ju = min(i__2,*n);
	mm = m;
	if (ju < kp1) {
	    goto L90;
	}
	i__2 = ju;
	for (j = kp1; j <= i__2; ++j) {
	    --l;
	    --mm;
	    t = abd[l + j * abd_dim1];
	    if (l == mm) {
		goto L70;
	    }
	    abd[l + j * abd_dim1] = abd[mm + j * abd_dim1];
	    abd[mm + j * abd_dim1] = t;
L70:
	    saxpy_(&lm, &t, &abd[m + 1 + k * abd_dim1], &c__1, &abd[mm + 1 + 
		    j * abd_dim1], &c__1);
/* L80: */
	}
L90:
	goto L110;
L100:
	*info = k;
L110:
/* L120: */
	;
    }
L130:
    ipvt[*n] = *n;
    if (abd[m + *n * abd_dim1] == 0.f) {
	*info = *n;
    }
    return 0;
} /* sgbfa_ */

#ifdef __cplusplus
}
#endif

