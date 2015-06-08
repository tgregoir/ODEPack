/* sstode.f -- translated by f2c (version 20100827).
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
    real conit, crate, el[13], elco[156]	/* was [13][12] */, hold, 
	    rmax, tesco[36]	/* was [3][12] */, ccmax, el0, h__, hmin, 
	    hmxi, hu, rc, tn, uround;
    integer iownd[6], ialth, ipup, lmax, meo, nqnyh, nslp, icf, ierpj, iersl, 
	    jcur, jstart, kflag, l, lyh, lewt, lacor, lsavf, lwm, liwm, meth, 
	    miter, maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu;
} sls001_;

#define sls001_1 sls001_

/* DECK SSTODE */
/* Subroutine */ int sstode_(integer *neq, real *y, real *yh, integer *nyh, 
	real *yh1, real *ewt, real *savf, real *acor, real *wm, integer *iwm, 
	S_fp f, U_fp jac, S_fp pjac, S_fp slvs)
{
    /* System generated locals */
    integer yh_dim1, yh_offset, i__1, i__2;
    real r__1, r__2, r__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j, m;
    static real r__;
    static integer i1, jb;
    static real rh, del, ddn;
    static integer ncf;
    static real dsm, dup, dcon, delp, rhdn, exdn;
    static integer iret;
    static real told, rhsm;
    static integer newq;
    static real exsm, rhup, exup;
    static integer iredo;
    extern /* Subroutine */ int scfode_(integer *, real *, real *);
    extern doublereal svnorm_(integer *, real *, real *);

/* ***BEGIN PROLOGUE  SSTODE */
/* ***SUBSIDIARY */
/* ***PURPOSE  Performs one step of an ODEPACK integration. */
/* ***TYPE      SINGLE PRECISION (SSTODE-S, DSTODE-D) */
/* ***AUTHOR  Hindmarsh, Alan C., (LLNL) */
/* ***DESCRIPTION */

/*  SSTODE performs one step of the integration of an initial value */
/*  problem for a system of ordinary differential equations. */
/*  Note:  SSTODE is independent of the value of the iteration method */
/*  indicator MITER, when this is .ne. 0, and hence is independent */
/*  of the type of chord method used, or the Jacobian structure. */
/*  Communication with SSTODE is done with the following variables: */

/*  NEQ    = integer array containing problem size in NEQ(1), and */
/*           passed as the NEQ argument in all calls to F and JAC. */
/*  Y      = an array of length .ge. N used as the Y argument in */
/*           all calls to F and JAC. */
/*  YH     = an NYH by LMAX array containing the dependent variables */
/*           and their approximate scaled derivatives, where */
/*           LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate */
/*           j-th derivative of y(i), scaled by h**j/factorial(j) */
/*           (j = 0,1,...,NQ).  on entry for the first step, the first */
/*           two columns of YH must be set from the initial values. */
/*  NYH    = a constant integer .ge. N, the first dimension of YH. */
/*  YH1    = a one-dimensional array occupying the same space as YH. */
/*  EWT    = an array of length N containing multiplicative weights */
/*           for local error measurements.  Local errors in Y(i) are */
/*           compared to 1.0/EWT(i) in various error tests. */
/*  SAVF   = an array of working storage, of length N. */
/*           Also used for input of YH(*,MAXORD+2) when JSTART = -1 */
/*           and MAXORD .lt. the current order NQ. */
/*  ACOR   = a work array of length N, used for the accumulated */
/*           corrections.  On a successful return, ACOR(i) contains */
/*           the estimated one-step local error in Y(i). */
/*  WM,IWM = real and integer work arrays associated with matrix */
/*           operations in chord iteration (MITER .ne. 0). */
/*  PJAC   = name of routine to evaluate and preprocess Jacobian matrix */
/*           and P = I - h*el0*JAC, if a chord method is being used. */
/*  SLVS   = name of routine to solve linear system in chord iteration. */
/*  CCMAX  = maximum relative change in h*el0 before PJAC is called. */
/*  H      = the step size to be attempted on the next step. */
/*           H is altered by the error control algorithm during the */
/*           problem.  H can be either positive or negative, but its */
/*           sign must remain constant throughout the problem. */
/*  HMIN   = the minimum absolute value of the step size h to be used. */
/*  HMXI   = inverse of the maximum absolute value of h to be used. */
/*           HMXI = 0.0 is allowed and corresponds to an infinite hmax. */
/*           HMIN and HMXI may be changed at any time, but will not */
/*           take effect until the next change of h is considered. */
/*  TN     = the independent variable. TN is updated on each step taken. */
/*  JSTART = an integer used for input only, with the following */
/*           values and meanings: */
/*                0  perform the first step. */
/*            .gt.0  take a new step continuing from the last. */
/*               -1  take the next step with a new value of H, MAXORD, */
/*                     N, METH, MITER, and/or matrix parameters. */
/*               -2  take the next step with a new value of H, */
/*                     but with other inputs unchanged. */
/*           On return, JSTART is set to 1 to facilitate continuation. */
/*  KFLAG  = a completion code with the following meanings: */
/*                0  the step was succesful. */
/*               -1  the requested error could not be achieved. */
/*               -2  corrector convergence could not be achieved. */
/*               -3  fatal error in PJAC or SLVS. */
/*           A return with KFLAG = -1 or -2 means either */
/*           abs(H) = HMIN or 10 consecutive failures occurred. */
/*           On a return with KFLAG negative, the values of TN and */
/*           the YH array are as of the beginning of the last */
/*           step, and H is the last step size attempted. */
/*  MAXORD = the maximum order of integration method to be allowed. */
/*  MAXCOR = the maximum number of corrector iterations allowed. */
/*  MSBP   = maximum number of steps between PJAC calls (MITER .gt. 0). */
/*  MXNCF  = maximum number of convergence failures allowed. */
/*  METH/MITER = the method flags.  See description in driver. */
/*  N      = the number of first-order differential equations. */
/*  The values of CCMAX, H, HMIN, HMXI, TN, JSTART, KFLAG, MAXORD, */
/*  MAXCOR, MSBP, MXNCF, METH, MITER, and N are communicated via COMMON. */

/* ***SEE ALSO  SLSODE */
/* ***ROUTINES CALLED  SCFODE, SVNORM */
/* ***COMMON BLOCKS    SLS001 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791129  DATE WRITTEN */
/*   890501  Modified prologue to SLATEC/LDOC format.  (FNF) */
/*   890503  Minor cosmetic changes.  (FNF) */
/*   930809  Renamed to allow single/double precision versions. (ACH) */
/*   010413  Reduced size of Common block /SLS001/. (ACH) */
/*   031105  Restored 'own' variables to Common block /SLS001/, to */
/*           enable interrupt/restart feature. (ACH) */
/* ***END PROLOGUE  SSTODE */
/* **End */

/* ***FIRST EXECUTABLE STATEMENT  SSTODE */
    /* Parameter adjustments */
    --neq;
    --y;
    yh_dim1 = *nyh;
    yh_offset = 1 + yh_dim1;
    yh -= yh_offset;
    --yh1;
    --ewt;
    --savf;
    --acor;
    --wm;
    --iwm;

    /* Function Body */
    sls001_1.kflag = 0;
    told = sls001_1.tn;
    ncf = 0;
    sls001_1.ierpj = 0;
    sls001_1.iersl = 0;
    sls001_1.jcur = 0;
    sls001_1.icf = 0;
    delp = 0.f;
    if (sls001_1.jstart > 0) {
	goto L200;
    }
    if (sls001_1.jstart == -1) {
	goto L100;
    }
    if (sls001_1.jstart == -2) {
	goto L160;
    }
/* ----------------------------------------------------------------------- */
/* On the first call, the order is set to 1, and other variables are */
/* initialized.  RMAX is the maximum ratio by which H can be increased */
/* in a single step.  It is initially 1.E4 to compensate for the small */
/* initial H, but then is normally equal to 10.  If a failure */
/* occurs (in corrector convergence or error test), RMAX is set to 2 */
/* for the next increase. */
/* ----------------------------------------------------------------------- */
    sls001_1.lmax = sls001_1.maxord + 1;
    sls001_1.nq = 1;
    sls001_1.l = 2;
    sls001_1.ialth = 2;
    sls001_1.rmax = 1e4f;
    sls001_1.rc = 0.f;
    sls001_1.el0 = 1.f;
    sls001_1.crate = .7f;
    sls001_1.hold = sls001_1.h__;
    sls001_1.meo = sls001_1.meth;
    sls001_1.nslp = 0;
    sls001_1.ipup = sls001_1.miter;
    iret = 3;
    goto L140;
/* ----------------------------------------------------------------------- */
/* The following block handles preliminaries needed when JSTART = -1. */
/* IPUP is set to MITER to force a matrix update. */
/* If an order increase is about to be considered (IALTH = 1), */
/* IALTH is reset to 2 to postpone consideration one more step. */
/* If the caller has changed METH, SCFODE is called to reset */
/* the coefficients of the method. */
/* If the caller has changed MAXORD to a value less than the current */
/* order NQ, NQ is reduced to MAXORD, and a new H chosen accordingly. */
/* If H is to be changed, YH must be rescaled. */
/* If H or METH is being changed, IALTH is reset to L = NQ + 1 */
/* to prevent further changes in H for that many steps. */
/* ----------------------------------------------------------------------- */
L100:
    sls001_1.ipup = sls001_1.miter;
    sls001_1.lmax = sls001_1.maxord + 1;
    if (sls001_1.ialth == 1) {
	sls001_1.ialth = 2;
    }
    if (sls001_1.meth == sls001_1.meo) {
	goto L110;
    }
    scfode_(&sls001_1.meth, sls001_1.elco, sls001_1.tesco);
    sls001_1.meo = sls001_1.meth;
    if (sls001_1.nq > sls001_1.maxord) {
	goto L120;
    }
    sls001_1.ialth = sls001_1.l;
    iret = 1;
    goto L150;
L110:
    if (sls001_1.nq <= sls001_1.maxord) {
	goto L160;
    }
L120:
    sls001_1.nq = sls001_1.maxord;
    sls001_1.l = sls001_1.lmax;
    i__1 = sls001_1.l;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L125: */
	sls001_1.el[i__ - 1] = sls001_1.elco[i__ + sls001_1.nq * 13 - 14];
    }
    sls001_1.nqnyh = sls001_1.nq * *nyh;
    sls001_1.rc = sls001_1.rc * sls001_1.el[0] / sls001_1.el0;
    sls001_1.el0 = sls001_1.el[0];
    sls001_1.conit = .5f / (sls001_1.nq + 2);
    ddn = svnorm_(&sls001_1.n, &savf[1], &ewt[1]) / sls001_1.tesco[sls001_1.l 
	    * 3 - 3];
    exdn = 1.f / sls001_1.l;
    d__1 = (doublereal) ddn;
    d__2 = (doublereal) exdn;
    rhdn = 1.f / (pow_dd(&d__1, &d__2) * 1.3f + 1.3e-6f);
    rh = dmin(rhdn,1.f);
    iredo = 3;
    if (sls001_1.h__ == sls001_1.hold) {
	goto L170;
    }
/* Computing MIN */
    r__2 = rh, r__3 = (r__1 = sls001_1.h__ / sls001_1.hold, dabs(r__1));
    rh = dmin(r__2,r__3);
    sls001_1.h__ = sls001_1.hold;
    goto L175;
/* ----------------------------------------------------------------------- */
/* SCFODE is called to get all the integration coefficients for the */
/* current METH.  Then the EL vector and related constants are reset */
/* whenever the order NQ is changed, or at the start of the problem. */
/* ----------------------------------------------------------------------- */
L140:
    scfode_(&sls001_1.meth, sls001_1.elco, sls001_1.tesco);
L150:
    i__1 = sls001_1.l;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L155: */
	sls001_1.el[i__ - 1] = sls001_1.elco[i__ + sls001_1.nq * 13 - 14];
    }
    sls001_1.nqnyh = sls001_1.nq * *nyh;
    sls001_1.rc = sls001_1.rc * sls001_1.el[0] / sls001_1.el0;
    sls001_1.el0 = sls001_1.el[0];
    sls001_1.conit = .5f / (sls001_1.nq + 2);
    switch (iret) {
	case 1:  goto L160;
	case 2:  goto L170;
	case 3:  goto L200;
    }
/* ----------------------------------------------------------------------- */
/* If H is being changed, the H ratio RH is checked against */
/* RMAX, HMIN, and HMXI, and the YH array rescaled.  IALTH is set to */
/* L = NQ + 1 to prevent a change of H for that many steps, unless */
/* forced by a convergence or error test failure. */
/* ----------------------------------------------------------------------- */
L160:
    if (sls001_1.h__ == sls001_1.hold) {
	goto L200;
    }
    rh = sls001_1.h__ / sls001_1.hold;
    sls001_1.h__ = sls001_1.hold;
    iredo = 3;
    goto L175;
L170:
/* Computing MAX */
    r__1 = rh, r__2 = sls001_1.hmin / dabs(sls001_1.h__);
    rh = dmax(r__1,r__2);
L175:
    rh = dmin(rh,sls001_1.rmax);
/* Computing MAX */
    r__1 = 1.f, r__2 = dabs(sls001_1.h__) * sls001_1.hmxi * rh;
    rh /= dmax(r__1,r__2);
    r__ = 1.f;
    i__1 = sls001_1.l;
    for (j = 2; j <= i__1; ++j) {
	r__ *= rh;
	i__2 = sls001_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L180: */
	    yh[i__ + j * yh_dim1] *= r__;
	}
    }
    sls001_1.h__ *= rh;
    sls001_1.rc *= rh;
    sls001_1.ialth = sls001_1.l;
    if (iredo == 0) {
	goto L690;
    }
/* ----------------------------------------------------------------------- */
/* This section computes the predicted values by effectively */
/* multiplying the YH array by the Pascal Triangle matrix. */
/* RC is the ratio of new to old values of the coefficient  H*EL(1). */
/* When RC differs from 1 by more than CCMAX, IPUP is set to MITER */
/* to force PJAC to be called, if a Jacobian is involved. */
/* In any case, PJAC is called at least every MSBP steps. */
/* ----------------------------------------------------------------------- */
L200:
    if ((r__1 = sls001_1.rc - 1.f, dabs(r__1)) > sls001_1.ccmax) {
	sls001_1.ipup = sls001_1.miter;
    }
    if (sls001_1.nst >= sls001_1.nslp + sls001_1.msbp) {
	sls001_1.ipup = sls001_1.miter;
    }
    sls001_1.tn += sls001_1.h__;
    i1 = sls001_1.nqnyh + 1;
    i__2 = sls001_1.nq;
    for (jb = 1; jb <= i__2; ++jb) {
	i1 -= *nyh;
/* dir$ ivdep */
	i__1 = sls001_1.nqnyh;
	for (i__ = i1; i__ <= i__1; ++i__) {
/* L210: */
	    yh1[i__] += yh1[i__ + *nyh];
	}
/* L215: */
    }
/* ----------------------------------------------------------------------- */
/* Up to MAXCOR corrector iterations are taken.  A convergence test is */
/* made on the R.M.S. norm of each correction, weighted by the error */
/* weight vector EWT.  The sum of the corrections is accumulated in the */
/* vector ACOR(i).  The YH array is not altered in the corrector loop. */
/* ----------------------------------------------------------------------- */
L220:
    m = 0;
    i__2 = sls001_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L230: */
	y[i__] = yh[i__ + yh_dim1];
    }
    (*f)(&neq[1], &sls001_1.tn, &y[1], &savf[1]);
    ++sls001_1.nfe;
    if (sls001_1.ipup <= 0) {
	goto L250;
    }
/* ----------------------------------------------------------------------- */
/* If indicated, the matrix P = I - h*el(1)*J is reevaluated and */
/* preprocessed before starting the corrector iteration.  IPUP is set */
/* to 0 as an indicator that this has been done. */
/* ----------------------------------------------------------------------- */
    (*pjac)(&neq[1], &y[1], &yh[yh_offset], nyh, &ewt[1], &acor[1], &savf[1], 
	    &wm[1], &iwm[1], (S_fp)f, (U_fp)jac);
    sls001_1.ipup = 0;
    sls001_1.rc = 1.f;
    sls001_1.nslp = sls001_1.nst;
    sls001_1.crate = .7f;
    if (sls001_1.ierpj != 0) {
	goto L430;
    }
L250:
    i__2 = sls001_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L260: */
	acor[i__] = 0.f;
    }
L270:
    if (sls001_1.miter != 0) {
	goto L350;
    }
/* ----------------------------------------------------------------------- */
/* In the case of functional iteration, update Y directly from */
/* the result of the last function evaluation. */
/* ----------------------------------------------------------------------- */
    i__2 = sls001_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	savf[i__] = sls001_1.h__ * savf[i__] - yh[i__ + (yh_dim1 << 1)];
/* L290: */
	y[i__] = savf[i__] - acor[i__];
    }
    del = svnorm_(&sls001_1.n, &y[1], &ewt[1]);
    i__2 = sls001_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	y[i__] = yh[i__ + yh_dim1] + sls001_1.el[0] * savf[i__];
/* L300: */
	acor[i__] = savf[i__];
    }
    goto L400;
/* ----------------------------------------------------------------------- */
/* In the case of the chord method, compute the corrector error, */
/* and solve the linear system with that as right-hand side and */
/* P as coefficient matrix. */
/* ----------------------------------------------------------------------- */
L350:
    i__2 = sls001_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L360: */
	y[i__] = sls001_1.h__ * savf[i__] - (yh[i__ + (yh_dim1 << 1)] + acor[
		i__]);
    }
    (*slvs)(&wm[1], &iwm[1], &y[1], &savf[1]);
    if (sls001_1.iersl < 0) {
	goto L430;
    }
    if (sls001_1.iersl > 0) {
	goto L410;
    }
    del = svnorm_(&sls001_1.n, &y[1], &ewt[1]);
    i__2 = sls001_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	acor[i__] += y[i__];
/* L380: */
	y[i__] = yh[i__ + yh_dim1] + sls001_1.el[0] * acor[i__];
    }
/* ----------------------------------------------------------------------- */
/* Test for convergence.  If M.gt.0, an estimate of the convergence */
/* rate constant is stored in CRATE, and this is used in the test. */
/* ----------------------------------------------------------------------- */
L400:
    if (m != 0) {
/* Computing MAX */
	r__1 = sls001_1.crate * .2f, r__2 = del / delp;
	sls001_1.crate = dmax(r__1,r__2);
    }
/* Computing MIN */
    r__1 = 1.f, r__2 = sls001_1.crate * 1.5f;
    dcon = del * dmin(r__1,r__2) / (sls001_1.tesco[sls001_1.nq * 3 - 2] * 
	    sls001_1.conit);
    if (dcon <= 1.f) {
	goto L450;
    }
    ++m;
    if (m == sls001_1.maxcor) {
	goto L410;
    }
    if (m >= 2 && del > delp * 2.f) {
	goto L410;
    }
    delp = del;
    (*f)(&neq[1], &sls001_1.tn, &y[1], &savf[1]);
    ++sls001_1.nfe;
    goto L270;
/* ----------------------------------------------------------------------- */
/* The corrector iteration failed to converge. */
/* If MITER .ne. 0 and the Jacobian is out of date, PJAC is called for */
/* the next try.  Otherwise the YH array is retracted to its values */
/* before prediction, and H is reduced, if possible.  If H cannot be */
/* reduced or MXNCF failures have occurred, exit with KFLAG = -2. */
/* ----------------------------------------------------------------------- */
L410:
    if (sls001_1.miter == 0 || sls001_1.jcur == 1) {
	goto L430;
    }
    sls001_1.icf = 1;
    sls001_1.ipup = sls001_1.miter;
    goto L220;
L430:
    sls001_1.icf = 2;
    ++ncf;
    sls001_1.rmax = 2.f;
    sls001_1.tn = told;
    i1 = sls001_1.nqnyh + 1;
    i__2 = sls001_1.nq;
    for (jb = 1; jb <= i__2; ++jb) {
	i1 -= *nyh;
/* dir$ ivdep */
	i__1 = sls001_1.nqnyh;
	for (i__ = i1; i__ <= i__1; ++i__) {
/* L440: */
	    yh1[i__] -= yh1[i__ + *nyh];
	}
/* L445: */
    }
    if (sls001_1.ierpj < 0 || sls001_1.iersl < 0) {
	goto L680;
    }
    if (dabs(sls001_1.h__) <= sls001_1.hmin * 1.00001f) {
	goto L670;
    }
    if (ncf == sls001_1.mxncf) {
	goto L670;
    }
    rh = .25f;
    sls001_1.ipup = sls001_1.miter;
    iredo = 1;
    goto L170;
/* ----------------------------------------------------------------------- */
/* The corrector has converged.  JCUR is set to 0 */
/* to signal that the Jacobian involved may need updating later. */
/* The local error test is made and control passes to statement 500 */
/* if it fails. */
/* ----------------------------------------------------------------------- */
L450:
    sls001_1.jcur = 0;
    if (m == 0) {
	dsm = del / sls001_1.tesco[sls001_1.nq * 3 - 2];
    }
    if (m > 0) {
	dsm = svnorm_(&sls001_1.n, &acor[1], &ewt[1]) / sls001_1.tesco[
		sls001_1.nq * 3 - 2];
    }
    if (dsm > 1.f) {
	goto L500;
    }
/* ----------------------------------------------------------------------- */
/* After a successful step, update the YH array. */
/* Consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1. */
/* If IALTH is then 1 and NQ .lt. MAXORD, then ACOR is saved for */
/* use in a possible order increase on the next step. */
/* If a change in H is considered, an increase or decrease in order */
/* by one is considered also.  A change in H is made only if it is by a */
/* factor of at least 1.1.  If not, IALTH is set to 3 to prevent */
/* testing for that many steps. */
/* ----------------------------------------------------------------------- */
    sls001_1.kflag = 0;
    iredo = 0;
    ++sls001_1.nst;
    sls001_1.hu = sls001_1.h__;
    sls001_1.nqu = sls001_1.nq;
    i__2 = sls001_1.l;
    for (j = 1; j <= i__2; ++j) {
	i__1 = sls001_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L470: */
	    yh[i__ + j * yh_dim1] += sls001_1.el[j - 1] * acor[i__];
	}
    }
    --sls001_1.ialth;
    if (sls001_1.ialth == 0) {
	goto L520;
    }
    if (sls001_1.ialth > 1) {
	goto L700;
    }
    if (sls001_1.l == sls001_1.lmax) {
	goto L700;
    }
    i__1 = sls001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L490: */
	yh[i__ + sls001_1.lmax * yh_dim1] = acor[i__];
    }
    goto L700;
/* ----------------------------------------------------------------------- */
/* The error test failed.  KFLAG keeps track of multiple failures. */
/* Restore TN and the YH array to their previous values, and prepare */
/* to try the step again.  Compute the optimum step size for this or */
/* one lower order.  After 2 or more failures, H is forced to decrease */
/* by a factor of 0.2 or less. */
/* ----------------------------------------------------------------------- */
L500:
    --sls001_1.kflag;
    sls001_1.tn = told;
    i1 = sls001_1.nqnyh + 1;
    i__1 = sls001_1.nq;
    for (jb = 1; jb <= i__1; ++jb) {
	i1 -= *nyh;
/* dir$ ivdep */
	i__2 = sls001_1.nqnyh;
	for (i__ = i1; i__ <= i__2; ++i__) {
/* L510: */
	    yh1[i__] -= yh1[i__ + *nyh];
	}
/* L515: */
    }
    sls001_1.rmax = 2.f;
    if (dabs(sls001_1.h__) <= sls001_1.hmin * 1.00001f) {
	goto L660;
    }
    if (sls001_1.kflag <= -3) {
	goto L640;
    }
    iredo = 2;
    rhup = 0.f;
    goto L540;
/* ----------------------------------------------------------------------- */
/* Regardless of the success or failure of the step, factors */
/* RHDN, RHSM, and RHUP are computed, by which H could be multiplied */
/* at order NQ - 1, order NQ, or order NQ + 1, respectively. */
/* In the case of failure, RHUP = 0.0 to avoid an order increase. */
/* The largest of these is determined and the new order chosen */
/* accordingly.  If the order is to be increased, we compute one */
/* additional scaled derivative. */
/* ----------------------------------------------------------------------- */
L520:
    rhup = 0.f;
    if (sls001_1.l == sls001_1.lmax) {
	goto L540;
    }
    i__1 = sls001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L530: */
	savf[i__] = acor[i__] - yh[i__ + sls001_1.lmax * yh_dim1];
    }
    dup = svnorm_(&sls001_1.n, &savf[1], &ewt[1]) / sls001_1.tesco[
	    sls001_1.nq * 3 - 1];
    exup = 1.f / (sls001_1.l + 1);
    d__1 = (doublereal) dup;
    d__2 = (doublereal) exup;
    rhup = 1.f / (pow_dd(&d__1, &d__2) * 1.4f + 1.4e-6f);
L540:
    exsm = 1.f / sls001_1.l;
    d__1 = (doublereal) dsm;
    d__2 = (doublereal) exsm;
    rhsm = 1.f / (pow_dd(&d__1, &d__2) * 1.2f + 1.2e-6f);
    rhdn = 0.f;
    if (sls001_1.nq == 1) {
	goto L560;
    }
    ddn = svnorm_(&sls001_1.n, &yh[sls001_1.l * yh_dim1 + 1], &ewt[1]) / 
	    sls001_1.tesco[sls001_1.nq * 3 - 3];
    exdn = 1.f / sls001_1.nq;
    d__1 = (doublereal) ddn;
    d__2 = (doublereal) exdn;
    rhdn = 1.f / (pow_dd(&d__1, &d__2) * 1.3f + 1.3e-6f);
L560:
    if (rhsm >= rhup) {
	goto L570;
    }
    if (rhup > rhdn) {
	goto L590;
    }
    goto L580;
L570:
    if (rhsm < rhdn) {
	goto L580;
    }
    newq = sls001_1.nq;
    rh = rhsm;
    goto L620;
L580:
    newq = sls001_1.nq - 1;
    rh = rhdn;
    if (sls001_1.kflag < 0 && rh > 1.f) {
	rh = 1.f;
    }
    goto L620;
L590:
    newq = sls001_1.l;
    rh = rhup;
    if (rh < 1.1f) {
	goto L610;
    }
    r__ = sls001_1.el[sls001_1.l - 1] / sls001_1.l;
    i__1 = sls001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L600: */
	yh[i__ + (newq + 1) * yh_dim1] = acor[i__] * r__;
    }
    goto L630;
L610:
    sls001_1.ialth = 3;
    goto L700;
L620:
    if (sls001_1.kflag == 0 && rh < 1.1f) {
	goto L610;
    }
    if (sls001_1.kflag <= -2) {
	rh = dmin(rh,.2f);
    }
/* ----------------------------------------------------------------------- */
/* If there is a change of order, reset NQ, l, and the coefficients. */
/* In any case H is reset according to RH and the YH array is rescaled. */
/* Then exit from 690 if the step was OK, or redo the step otherwise. */
/* ----------------------------------------------------------------------- */
    if (newq == sls001_1.nq) {
	goto L170;
    }
L630:
    sls001_1.nq = newq;
    sls001_1.l = sls001_1.nq + 1;
    iret = 2;
    goto L150;
/* ----------------------------------------------------------------------- */
/* Control reaches this section if 3 or more failures have occured. */
/* If 10 failures have occurred, exit with KFLAG = -1. */
/* It is assumed that the derivatives that have accumulated in the */
/* YH array have errors of the wrong order.  Hence the first */
/* derivative is recomputed, and the order is set to 1.  Then */
/* H is reduced by a factor of 10, and the step is retried, */
/* until it succeeds or H reaches HMIN. */
/* ----------------------------------------------------------------------- */
L640:
    if (sls001_1.kflag == -10) {
	goto L660;
    }
    rh = .1f;
/* Computing MAX */
    r__1 = sls001_1.hmin / dabs(sls001_1.h__);
    rh = dmax(r__1,rh);
    sls001_1.h__ *= rh;
    i__1 = sls001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L645: */
	y[i__] = yh[i__ + yh_dim1];
    }
    (*f)(&neq[1], &sls001_1.tn, &y[1], &savf[1]);
    ++sls001_1.nfe;
    i__1 = sls001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L650: */
	yh[i__ + (yh_dim1 << 1)] = sls001_1.h__ * savf[i__];
    }
    sls001_1.ipup = sls001_1.miter;
    sls001_1.ialth = 5;
    if (sls001_1.nq == 1) {
	goto L200;
    }
    sls001_1.nq = 1;
    sls001_1.l = 2;
    iret = 3;
    goto L150;
/* ----------------------------------------------------------------------- */
/* All returns are made through this section.  H is saved in HOLD */
/* to allow the caller to change H on the next step. */
/* ----------------------------------------------------------------------- */
L660:
    sls001_1.kflag = -1;
    goto L720;
L670:
    sls001_1.kflag = -2;
    goto L720;
L680:
    sls001_1.kflag = -3;
    goto L720;
L690:
    sls001_1.rmax = 10.f;
L700:
    r__ = 1.f / sls001_1.tesco[sls001_1.nqu * 3 - 2];
    i__1 = sls001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L710: */
	acor[i__] *= r__;
    }
L720:
    sls001_1.hold = sls001_1.h__;
    sls001_1.jstart = 1;
    return 0;
/* ----------------------- END OF SUBROUTINE SSTODE ---------------------- */
} /* sstode_ */

