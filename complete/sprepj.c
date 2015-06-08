/* sprepj.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    real rowns[209], ccmax, el0, h__, hmin, hmxi, hu, rc, tn, uround;
    integer iownd[6], iowns[6], icf, ierpj, iersl, jcur, jstart, kflag, l, 
	    lyh, lewt, lacor, lsavf, lwm, liwm, meth, miter, maxord, maxcor, 
	    msbp, mxncf, n, nq, nst, nfe, nje, nqu;
} sls001_;

#define sls001_1 sls001_

/* Table of constant values */

static integer c__0 = 0;

/* DECK SPREPJ */
/* Subroutine */ int sprepj_(integer *neq, real *y, real *yh, integer *nyh, 
	real *ewt, real *ftem, real *savf, real *wm, integer *iwm, S_fp f, 
	S_fp jac)
{
    /* System generated locals */
    integer yh_dim1, yh_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2;

    /* Local variables */
    static integer i__, j;
    static real r__;
    static integer i1, i2, j1;
    static real r0, di;
    static integer ii, jj, ml, mu;
    static real yi, yj, hl0;
    static integer ml3, np1;
    static real fac;
    static integer mba, ier;
    static real con, yjj;
    static integer meb1, lenp;
    static real srur;
    static integer mband;
    extern /* Subroutine */ int sgbfa_(real *, integer *, integer *, integer *
	    , integer *, integer *, integer *), sgefa_(real *, integer *, 
	    integer *, integer *, integer *);
    static integer meband;
    extern doublereal svnorm_(integer *, real *, real *);

/* ***BEGIN PROLOGUE  SPREPJ */
/* ***SUBSIDIARY */
/* ***PURPOSE  Compute and process Newton iteration matrix. */
/* ***TYPE      SINGLE PRECISION (SPREPJ-S, DPREPJ-D) */
/* ***AUTHOR  Hindmarsh, Alan C., (LLNL) */
/* ***DESCRIPTION */

/*  SPREPJ is called by SSTODE to compute and process the matrix */
/*  P = I - h*el(1)*J , where J is an approximation to the Jacobian. */
/*  Here J is computed by the user-supplied routine JAC if */
/*  MITER = 1 or 4, or by finite differencing if MITER = 2, 3, or 5. */
/*  If MITER = 3, a diagonal approximation to J is used. */
/*  J is stored in WM and replaced by P.  If MITER .ne. 3, P is then */
/*  subjected to LU decomposition in preparation for later solution */
/*  of linear systems with P as coefficient matrix.  This is done */
/*  by SGEFA if MITER = 1 or 2, and by SGBFA if MITER = 4 or 5. */

/*  In addition to variables described in SSTODE and SLSODE prologues, */
/*  communication with SPREPJ uses the following: */
/*  Y     = array containing predicted values on entry. */
/*  FTEM  = work array of length N (ACOR in SSTODE). */
/*  SAVF  = array containing f evaluated at predicted y. */
/*  WM    = real work space for matrices.  On output it contains the */
/*          inverse diagonal matrix if MITER = 3 and the LU decomposition */
/*          of P if MITER is 1, 2 , 4, or 5. */
/*          Storage of matrix elements starts at WM(3). */
/*          WM also contains the following matrix-related data: */
/*          WM(1) = SQRT(UROUND), used in numerical Jacobian increments. */
/*          WM(2) = H*EL0, saved for later use if MITER = 3. */
/*  IWM   = integer work space containing pivot information, starting at */
/*          IWM(21), if MITER is 1, 2, 4, or 5.  IWM also contains band */
/*          parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5. */
/*  EL0   = EL(1) (input). */
/*  IERPJ = output error flag,  = 0 if no trouble, .gt. 0 if */
/*          P matrix found to be singular. */
/*  JCUR  = output flag = 1 to indicate that the Jacobian matrix */
/*          (or approximation) is now current. */
/*  This routine also uses the COMMON variables EL0, H, TN, UROUND, */
/*  MITER, N, NFE, and NJE. */

/* ***SEE ALSO  SLSODE */
/* ***ROUTINES CALLED  SGBFA, SGEFA, SVNORM */
/* ***COMMON BLOCKS    SLS001 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791129  DATE WRITTEN */
/*   890501  Modified prologue to SLATEC/LDOC format.  (FNF) */
/*   890504  Minor cosmetic changes.  (FNF) */
/*   930809  Renamed to allow single/double precision versions. (ACH) */
/*   010412  Reduced size of Common block /SLS001/. (ACH) */
/*   031105  Restored 'own' variables to Common block /SLS001/, to */
/*           enable interrupt/restart feature. (ACH) */
/* ***END PROLOGUE  SPREPJ */
/* **End */

/* ***FIRST EXECUTABLE STATEMENT  SPREPJ */
    /* Parameter adjustments */
    --neq;
    --y;
    yh_dim1 = *nyh;
    yh_offset = 1 + yh_dim1;
    yh -= yh_offset;
    --ewt;
    --ftem;
    --savf;
    --wm;
    --iwm;

    /* Function Body */
    ++sls001_1.nje;
    sls001_1.ierpj = 0;
    sls001_1.jcur = 1;
    hl0 = sls001_1.h__ * sls001_1.el0;
    switch (sls001_1.miter) {
	case 1:  goto L100;
	case 2:  goto L200;
	case 3:  goto L300;
	case 4:  goto L400;
	case 5:  goto L500;
    }
/* If MITER = 1, call JAC and multiply by scalar. ----------------------- */
L100:
    lenp = sls001_1.n * sls001_1.n;
    i__1 = lenp;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L110: */
	wm[i__ + 2] = 0.f;
    }
    (*jac)(&neq[1], &sls001_1.tn, &y[1], &c__0, &c__0, &wm[3], &sls001_1.n);
    con = -hl0;
    i__1 = lenp;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L120: */
	wm[i__ + 2] *= con;
    }
    goto L240;
/* If MITER = 2, make N calls to F to approximate J. -------------------- */
L200:
    fac = svnorm_(&sls001_1.n, &savf[1], &ewt[1]);
    r0 = dabs(sls001_1.h__) * 1e3f * sls001_1.uround * sls001_1.n * fac;
    if (r0 == 0.f) {
	r0 = 1.f;
    }
    srur = wm[1];
    j1 = 2;
    i__1 = sls001_1.n;
    for (j = 1; j <= i__1; ++j) {
	yj = y[j];
/* Computing MAX */
	r__1 = srur * dabs(yj), r__2 = r0 / ewt[j];
	r__ = dmax(r__1,r__2);
	y[j] += r__;
	fac = -hl0 / r__;
	(*f)(&neq[1], &sls001_1.tn, &y[1], &ftem[1]);
	i__2 = sls001_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L220: */
	    wm[i__ + j1] = (ftem[i__] - savf[i__]) * fac;
	}
	y[j] = yj;
	j1 += sls001_1.n;
/* L230: */
    }
    sls001_1.nfe += sls001_1.n;
/* Add identity matrix. ------------------------------------------------- */
L240:
    j = 3;
    np1 = sls001_1.n + 1;
    i__1 = sls001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wm[j] += 1.f;
/* L250: */
	j += np1;
    }
/* Do LU decomposition on P. -------------------------------------------- */
    sgefa_(&wm[3], &sls001_1.n, &sls001_1.n, &iwm[21], &ier);
    if (ier != 0) {
	sls001_1.ierpj = 1;
    }
    return 0;
/* If MITER = 3, construct a diagonal approximation to J and P. --------- */
L300:
    wm[2] = hl0;
    r__ = sls001_1.el0 * .1f;
    i__1 = sls001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L310: */
	y[i__] += r__ * (sls001_1.h__ * savf[i__] - yh[i__ + (yh_dim1 << 1)]);
    }
    (*f)(&neq[1], &sls001_1.tn, &y[1], &wm[3]);
    ++sls001_1.nfe;
    i__1 = sls001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r0 = sls001_1.h__ * savf[i__] - yh[i__ + (yh_dim1 << 1)];
	di = r0 * .1f - sls001_1.h__ * (wm[i__ + 2] - savf[i__]);
	wm[i__ + 2] = 1.f;
	if (dabs(r0) < sls001_1.uround / ewt[i__]) {
	    goto L320;
	}
	if (dabs(di) == 0.f) {
	    goto L330;
	}
	wm[i__ + 2] = r0 * .1f / di;
L320:
	;
    }
    return 0;
L330:
    sls001_1.ierpj = 1;
    return 0;
/* If MITER = 4, call JAC and multiply by scalar. ----------------------- */
L400:
    ml = iwm[1];
    mu = iwm[2];
    ml3 = ml + 3;
    mband = ml + mu + 1;
    meband = mband + ml;
    lenp = meband * sls001_1.n;
    i__1 = lenp;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L410: */
	wm[i__ + 2] = 0.f;
    }
    (*jac)(&neq[1], &sls001_1.tn, &y[1], &ml, &mu, &wm[ml3], &meband);
    con = -hl0;
    i__1 = lenp;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L420: */
	wm[i__ + 2] *= con;
    }
    goto L570;
/* If MITER = 5, make MBAND calls to F to approximate J. ---------------- */
L500:
    ml = iwm[1];
    mu = iwm[2];
    mband = ml + mu + 1;
    mba = min(mband,sls001_1.n);
    meband = mband + ml;
    meb1 = meband - 1;
    srur = wm[1];
    fac = svnorm_(&sls001_1.n, &savf[1], &ewt[1]);
    r0 = dabs(sls001_1.h__) * 1e3f * sls001_1.uround * sls001_1.n * fac;
    if (r0 == 0.f) {
	r0 = 1.f;
    }
    i__1 = mba;
    for (j = 1; j <= i__1; ++j) {
	i__2 = sls001_1.n;
	i__3 = mband;
	for (i__ = j; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__3) {
	    yi = y[i__];
/* Computing MAX */
	    r__1 = srur * dabs(yi), r__2 = r0 / ewt[i__];
	    r__ = dmax(r__1,r__2);
/* L530: */
	    y[i__] += r__;
	}
	(*f)(&neq[1], &sls001_1.tn, &y[1], &ftem[1]);
	i__3 = sls001_1.n;
	i__2 = mband;
	for (jj = j; i__2 < 0 ? jj >= i__3 : jj <= i__3; jj += i__2) {
	    y[jj] = yh[jj + yh_dim1];
	    yjj = y[jj];
/* Computing MAX */
	    r__1 = srur * dabs(yjj), r__2 = r0 / ewt[jj];
	    r__ = dmax(r__1,r__2);
	    fac = -hl0 / r__;
/* Computing MAX */
	    i__4 = jj - mu;
	    i1 = max(i__4,1);
/* Computing MIN */
	    i__4 = jj + ml;
	    i2 = min(i__4,sls001_1.n);
	    ii = jj * meb1 - ml + 2;
	    i__4 = i2;
	    for (i__ = i1; i__ <= i__4; ++i__) {
/* L540: */
		wm[ii + i__] = (ftem[i__] - savf[i__]) * fac;
	    }
/* L550: */
	}
/* L560: */
    }
    sls001_1.nfe += mba;
/* Add identity matrix. ------------------------------------------------- */
L570:
    ii = mband + 2;
    i__1 = sls001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wm[ii] += 1.f;
/* L580: */
	ii += meband;
    }
/* Do LU decomposition of P. -------------------------------------------- */
    sgbfa_(&wm[3], &meband, &sls001_1.n, &ml, &mu, &iwm[21], &ier);
    if (ier != 0) {
	sls001_1.ierpj = 1;
    }
    return 0;
/* ----------------------- END OF SUBROUTINE SPREPJ ---------------------- */
} /* sprepj_ */

