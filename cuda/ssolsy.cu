#ifdef __cplusplus
extern "C" {
#endif

#include "f2c.h"

/* DECK SSOLSY */
/* Subroutine */
__device__ int ssolsy_(real *wm, integer *iwm, real *x, real *tem)
{
    /* Common Block Declarations */

    struct {
        real rowns[209], ccmax, el0, h__, hmin, hmxi, hu, rc, tn, uround;
        integer iownd[6], iowns[6], icf, ierpj, iersl, jcur, jstart, kflag, l, 
	        lyh, lewt, lacor, lsavf, lwm, liwm, meth, miter, maxord, maxcor, 
	        msbp, mxncf, n, nq, nst, nfe, nje, nqu;
    } sls001_;

    #define sls001_1 sls001_

    /* Table of constant values */

    integer c__0 = 0;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__;
    real r__, di;
    integer ml, mu;
    real hl0, phl0;
    extern /* Subroutine */ __device__ int sgbsl_(real *, integer *, integer *, integer *
	    , integer *, integer *, real *, integer *), sgesl_(real *, 
	    integer *, integer *, integer *, real *, integer *);
    integer meband;

/* ***BEGIN PROLOGUE  SSOLSY */
/* ***SUBSIDIARY */
/* ***PURPOSE  ODEPACK linear system solver. */
/* ***TYPE      SINGLE PRECISION (SSOLSY-S, DSOLSY-D) */
/* ***AUTHOR  Hindmarsh, Alan C., (LLNL) */
/* ***DESCRIPTION */

/*  This routine manages the solution of the linear system arising from */
/*  a chord iteration.  It is called if MITER .ne. 0. */
/*  If MITER is 1 or 2, it calls SGESL to accomplish this. */
/*  If MITER = 3 it updates the coefficient h*EL0 in the diagonal */
/*  matrix, and then computes the solution. */
/*  If MITER is 4 or 5, it calls SGBSL. */
/*  Communication with SSOLSY uses the following variables: */
/*  WM    = real work space containing the inverse diagonal matrix if */
/*          MITER = 3 and the LU decomposition of the matrix otherwise. */
/*          Storage of matrix elements starts at WM(3). */
/*          WM also contains the following matrix-related data: */
/*          WM(1) = SQRT(UROUND) (not used here), */
/*          WM(2) = HL0, the previous value of h*EL0, used if MITER = 3. */
/*  IWM   = integer work space containing pivot information, starting at */
/*          IWM(21), if MITER is 1, 2, 4, or 5.  IWM also contains band */
/*          parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5. */
/*  X     = the right-hand side vector on input, and the solution vector */
/*          on output, of length N. */
/*  TEM   = vector of work space of length N, not used in this version. */
/*  IERSL = output flag (in COMMON).  IERSL = 0 if no trouble occurred. */
/*          IERSL = 1 if a singular matrix arose with MITER = 3. */
/*  This routine also uses the COMMON variables EL0, H, MITER, and N. */

/* ***SEE ALSO  SLSODE */
/* ***ROUTINES CALLED  SGBSL, SGESL */
/* ***COMMON BLOCKS    SLS001 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791129  DATE WRITTEN */
/*   890501  Modified prologue to SLATEC/LDOC format.  (FNF) */
/*   890503  Minor cosmetic changes.  (FNF) */
/*   930809  Renamed to allow single/double precision versions. (ACH) */
/*   010412  Reduced size of Common block /SLS001/. (ACH) */
/*   031105  Restored 'own' variables to Common block /SLS001/, to */
/*           enable interrupt/restart feature. (ACH) */
/* ***END PROLOGUE  SSOLSY */
/* **End */

/* ***FIRST EXECUTABLE STATEMENT  SSOLSY */
    /* Parameter adjustments */
    --tem;
    --x;
    --iwm;
    --wm;

    /* Function Body */
    sls001_1.iersl = 0;
    switch (sls001_1.miter) {
	case 1:  goto L100;
	case 2:  goto L100;
	case 3:  goto L300;
	case 4:  goto L400;
	case 5:  goto L400;
    }
L100:
    sgesl_(&wm[3], &sls001_1.n, &sls001_1.n, &iwm[21], &x[1], &c__0);
    return 0;

L300:
    phl0 = wm[2];
    hl0 = sls001_1.h__ * sls001_1.el0;
    wm[2] = hl0;
    if (hl0 == phl0) {
	goto L330;
    }
    r__ = hl0 / phl0;
    i__1 = sls001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	di = 1.f - r__ * (1.f - 1.f / wm[i__ + 2]);
	if (dabs(di) == 0.f) {
	    goto L390;
	}
/* L320: */
	wm[i__ + 2] = 1.f / di;
    }
L330:
    i__1 = sls001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L340: */
	x[i__] = wm[i__ + 2] * x[i__];
    }
    return 0;
L390:
    sls001_1.iersl = 1;
    return 0;

L400:
    ml = iwm[1];
    mu = iwm[2];
    meband = (ml << 1) + mu + 1;
    sgbsl_(&wm[3], &meband, &sls001_1.n, &ml, &mu, &iwm[21], &x[1], &c__0);
    return 0;
/* ----------------------- END OF SUBROUTINE SSOLSY ---------------------- */
} /* ssolsy_ */

#ifdef __cplusplus
}
#endif

