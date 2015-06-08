#ifdef __cplusplus
extern "C" {
#endif

#include "f2c.h"

/* DECK SVNORM */
__device__ doublereal svnorm_(integer *n, real *v, real *w)
{
    /* System generated locals */
    integer i__1;
    real ret_val, r__1;

    /* Builtin functions */
    //double sqrt(doublereal);

    /* Local variables */
    integer i__;
    real sum;

/* ***BEGIN PROLOGUE  SVNORM */
/* ***SUBSIDIARY */
/* ***PURPOSE  Weighted root-mean-square vector norm. */
/* ***TYPE      SINGLE PRECISION (SVNORM-S, DVNORM-D) */
/* ***AUTHOR  Hindmarsh, Alan C., (LLNL) */
/* ***DESCRIPTION */

/*  This function routine computes the weighted root-mean-square norm */
/*  of the vector of length N contained in the array V, with weights */
/*  contained in the array W of length N: */
/*    SVNORM = SQRT( (1/N) * SUM( V(i)*W(i) )**2 ) */

/* ***SEE ALSO  SLSODE */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791129  DATE WRITTEN */
/*   890501  Modified prologue to SLATEC/LDOC format.  (FNF) */
/*   890503  Minor cosmetic changes.  (FNF) */
/*   930809  Renamed to allow single/double precision versions. (ACH) */
/* ***END PROLOGUE  SVNORM */
/* **End */

/* ***FIRST EXECUTABLE STATEMENT  SVNORM */
    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    sum = 0.f;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
/* Computing 2nd power */
	r__1 = v[i__] * w[i__];
	sum += r__1 * r__1;
    }
    ret_val = sqrt(sum / *n);
    return ret_val;
/* ----------------------- END OF FUNCTION SVNORM ------------------------ */
} /* svnorm_ */

#ifdef __cplusplus
}
#endif

