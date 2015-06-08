#ifdef __cplusplus
extern "C" {
#endif

#include "f2c.h"

/* DECK SINTDY */
/* Subroutine */
__device__ int sintdy_(real *t, integer *k, real *yh, integer *nyh, 
	real *dky, integer *iflag)
{
    /* Common Block Declarations */

    struct {
        real rowns[209], ccmax, el0, h__, hmin, hmxi, hu, rc, tn, uround;
        integer iownd[6], iowns[6], icf, ierpj, iersl, jcur, jstart, kflag, l, 
	        lyh, lewt, lacor, lsavf, lwm, liwm, meth, miter, maxord, maxcor, 
	        msbp, mxncf, n, nq, nst, nfe, nje, nqu;
    } sls001_;

    #define sls001_1 sls001_

    /* System generated locals */
    integer yh_dim1, yh_offset, i__1, i__2;
    real r__1;

    /* Builtin functions */
    extern __device__ double r_sign(real *, real *), pow_ri(real *, integer *);

    /* Local variables */
    real c__;
    integer i__, j;
    real r__, s;
    integer ic, jb, jj;
    real tp;
    integer jb2, jj1, jp1;

/* ***BEGIN PROLOGUE  SINTDY */
/* ***SUBSIDIARY */
/* ***PURPOSE  Interpolate solution derivatives. */
/* ***TYPE      SINGLE PRECISION (SINTDY-S, DINTDY-D) */
/* ***AUTHOR  Hindmarsh, Alan C., (LLNL) */
/* ***DESCRIPTION */

/*  SINTDY computes interpolated values of the K-th derivative of the */
/*  dependent variable vector y, and stores it in DKY.  This routine */
/*  is called within the package with K = 0 and T = TOUT, but may */
/*  also be called by the user for any K up to the current order. */
/*  (See detailed instructions in the usage documentation.) */

/*  The computed values in DKY are gotten by interpolation using the */
/*  Nordsieck history array YH.  This array corresponds uniquely to a */
/*  vector-valued polynomial of degree NQCUR or less, and DKY is set */
/*  to the K-th derivative of this polynomial at T. */
/*  The formula for DKY is: */
/*               q */
/*   DKY(i)  =  sum  c(j,K) * (T - tn)**(j-K) * h**(-j) * YH(i,j+1) */
/*              j=K */
/*  where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, tn = TCUR, h = HCUR. */
/*  The quantities  nq = NQCUR, l = nq+1, N = NEQ, tn, and h are */
/*  communicated by COMMON.  The above sum is done in reverse order. */
/*  IFLAG is returned negative if either K or T is out of bounds. */

/* ***SEE ALSO  SLSODE */
/* ***ROUTINES CALLED  XERRWV */
/* ***COMMON BLOCKS    SLS001 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791129  DATE WRITTEN */
/*   890501  Modified prologue to SLATEC/LDOC format.  (FNF) */
/*   890503  Minor cosmetic changes.  (FNF) */
/*   930809  Renamed to allow single/double precision versions. (ACH) */
/*   010412  Reduced size of Common block /SLS001/. (ACH) */
/*   031105  Restored 'own' variables to Common block /SLS001/, to */
/*           enable interrupt/restart feature. (ACH) */
/*   050427  Corrected roundoff decrement in TP. (ACH) */
/* ***END PROLOGUE  SINTDY */
/* **End */

/* ***FIRST EXECUTABLE STATEMENT  SINTDY */
    /* Parameter adjustments */
    yh_dim1 = *nyh;
    yh_offset = 1 + yh_dim1;
    yh -= yh_offset;
    --dky;

    /* Function Body */
    *iflag = 0;
    if (*k < 0 || *k > sls001_1.nq) {
	goto L80;
    }
    r__1 = dabs(sls001_1.tn) + dabs(sls001_1.hu);
    tp = sls001_1.tn - sls001_1.hu - sls001_1.uround * 100.f * r_sign(&r__1, &
	    sls001_1.hu);
    if ((*t - tp) * (*t - sls001_1.tn) > 0.f) {
	goto L90;
    }

    s = (*t - sls001_1.tn) / sls001_1.h__;
    ic = 1;
    if (*k == 0) {
	goto L15;
    }
    jj1 = sls001_1.l - *k;
    i__1 = sls001_1.nq;
    for (jj = jj1; jj <= i__1; ++jj) {
/* L10: */
	ic *= jj;
    }
L15:
    c__ = (real) ic;
    i__1 = sls001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L20: */
	dky[i__] = c__ * yh[i__ + sls001_1.l * yh_dim1];
    }
    if (*k == sls001_1.nq) {
	goto L55;
    }
    jb2 = sls001_1.nq - *k;
    i__1 = jb2;
    for (jb = 1; jb <= i__1; ++jb) {
	j = sls001_1.nq - jb;
	jp1 = j + 1;
	ic = 1;
	if (*k == 0) {
	    goto L35;
	}
	jj1 = jp1 - *k;
	i__2 = j;
	for (jj = jj1; jj <= i__2; ++jj) {
/* L30: */
	    ic *= jj;
	}
L35:
	c__ = (real) ic;
	i__2 = sls001_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L40: */
	    dky[i__] = c__ * yh[i__ + jp1 * yh_dim1] + s * dky[i__];
	}
/* L50: */
    }
    if (*k == 0) {
	return 0;
    }
L55:
    i__1 = -(*k);
    r__ = pow_ri(&sls001_1.h__, &i__1);
    i__1 = sls001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L60: */
	dky[i__] = r__ * dky[i__];
    }
    return 0;

L80:
    *iflag = -1;
    return 0;
L90:
    *iflag = -2;
    return 0;
/* ----------------------- END OF SUBROUTINE SINTDY ---------------------- */
} /* sintdy_ */

#ifdef __cplusplus
}
#endif

