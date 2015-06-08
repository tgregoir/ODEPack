#ifdef __cplusplus
extern "C" {
#endif

#include "f2c.h"

/* DECK SGEFA */
/* Subroutine */
__device__ int sgefa_(real *a, integer *lda, integer *n, integer *ipvt, 
	integer *info)
{
    /* Table of constant values */

    integer c__1 = 1;

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    integer j, k, l;
    real t;
    integer kp1, nm1;
    extern /* Subroutine */ __device__ int sscal_(integer *, real *, real *, integer *), 
	    saxpy_(integer *, real *, real *, integer *, real *, integer *);
    extern __device__ integer isamax_(integer *, real *, integer *);

/* ***BEGIN PROLOGUE  SGEFA */
/* ***PURPOSE  Factor a matrix using Gaussian elimination. */
/* ***CATEGORY  D2A1 */
/* ***TYPE      SINGLE PRECISION (SGEFA-S, DGEFA-D, CGEFA-C) */
/* ***KEYWORDS  GENERAL MATRIX, LINEAR ALGEBRA, LINPACK, */
/*             MATRIX FACTORIZATION */
/* ***AUTHOR  Moler, C. B., (U. of New Mexico) */
/* ***DESCRIPTION */

/*     SGEFA factors a real matrix by Gaussian elimination. */

/*     SGEFA is usually called by SGECO, but it can be called */
/*     directly with a saving in time if  RCOND  is not needed. */
/*     (Time for SGECO) = (1 + 9/N)*(Time for SGEFA) . */

/*     On Entry */

/*        A       REAL(LDA, N) */
/*                the matrix to be factored. */

/*        LDA     INTEGER */
/*                the leading dimension of the array  A . */

/*        N       INTEGER */
/*                the order of the matrix  A . */

/*     On Return */

/*        A       an upper triangular matrix and the multipliers */
/*                which were used to obtain it. */
/*                The factorization can be written  A = L*U , where */
/*                L  is a product of permutation and unit lower */
/*                triangular matrices and  U  is upper triangular. */

/*        IPVT    INTEGER(N) */
/*                an integer vector of pivot indices. */

/*        INFO    INTEGER */
/*                = 0  normal value. */
/*                = K  if  U(K,K) .EQ. 0.0 .  This is not an error */
/*                     condition for this subroutine, but it does */
/*                     indicate that SGESL or SGEDI will divide by zero */
/*                     if called.  Use  RCOND  in SGECO for a reliable */
/*                     indication of singularity. */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  ISAMAX, SAXPY, SSCAL */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  SGEFA */


/*     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING */

/* ***FIRST EXECUTABLE STATEMENT  SGEFA */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipvt;

    /* Function Body */
    *info = 0;
    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L70;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;

/*        FIND L = PIVOT INDEX */

	i__2 = *n - k + 1;
	l = isamax_(&i__2, &a[k + k * a_dim1], &c__1) + k - 1;
	ipvt[k] = l;

/*        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED */

	if (a[l + k * a_dim1] == 0.f) {
	    goto L40;
	}

/*           INTERCHANGE IF NECESSARY */

	if (l == k) {
	    goto L10;
	}
	t = a[l + k * a_dim1];
	a[l + k * a_dim1] = a[k + k * a_dim1];
	a[k + k * a_dim1] = t;
L10:

/*           COMPUTE MULTIPLIERS */

	t = -1.f / a[k + k * a_dim1];
	i__2 = *n - k;
	sscal_(&i__2, &t, &a[k + 1 + k * a_dim1], &c__1);

/*           ROW ELIMINATION WITH COLUMN INDEXING */

	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    t = a[l + j * a_dim1];
	    if (l == k) {
		goto L20;
	    }
	    a[l + j * a_dim1] = a[k + j * a_dim1];
	    a[k + j * a_dim1] = t;
L20:
	    i__3 = *n - k;
	    saxpy_(&i__3, &t, &a[k + 1 + k * a_dim1], &c__1, &a[k + 1 + j * 
		    a_dim1], &c__1);
/* L30: */
	}
	goto L50;
L40:
	*info = k;
L50:
/* L60: */
	;
    }
L70:
    ipvt[*n] = *n;
    if (a[*n + *n * a_dim1] == 0.f) {
	*info = *n;
    }
    return 0;
} /* sgefa_ */

#ifdef __cplusplus
}
#endif

