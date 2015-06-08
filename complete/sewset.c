/* sewset.f -- translated by f2c (version 20100827).
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

/* DECK SEWSET */
/* Subroutine */ int sewset_(integer *n, integer *itol, real *rtol, real *
	atol, real *ycur, real *ewt)
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Local variables */
    static integer i__;

/* ***BEGIN PROLOGUE  SEWSET */
/* ***SUBSIDIARY */
/* ***PURPOSE  Set error weight vector. */
/* ***TYPE      SINGLE PRECISION (SEWSET-S, DEWSET-D) */
/* ***AUTHOR  Hindmarsh, Alan C., (LLNL) */
/* ***DESCRIPTION */

/*  This subroutine sets the error weight vector EWT according to */
/*      EWT(i) = RTOL(i)*ABS(YCUR(i)) + ATOL(i),  i = 1,...,N, */
/*  with the subscript on RTOL and/or ATOL possibly replaced by 1 above, */
/*  depending on the value of ITOL. */

/* ***SEE ALSO  SLSODE */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791129  DATE WRITTEN */
/*   890501  Modified prologue to SLATEC/LDOC format.  (FNF) */
/*   890503  Minor cosmetic changes.  (FNF) */
/*   930809  Renamed to allow single/double precision versions. (ACH) */
/* ***END PROLOGUE  SEWSET */
/* **End */

/* ***FIRST EXECUTABLE STATEMENT  SEWSET */
    /* Parameter adjustments */
    --ewt;
    --ycur;
    --rtol;
    --atol;

    /* Function Body */
    switch (*itol) {
	case 1:  goto L10;
	case 2:  goto L20;
	case 3:  goto L30;
	case 4:  goto L40;
    }
L10:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L15: */
	ewt[i__] = rtol[1] * (r__1 = ycur[i__], dabs(r__1)) + atol[1];
    }
    return 0;
L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L25: */
	ewt[i__] = rtol[1] * (r__1 = ycur[i__], dabs(r__1)) + atol[i__];
    }
    return 0;
L30:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L35: */
	ewt[i__] = rtol[i__] * (r__1 = ycur[i__], dabs(r__1)) + atol[1];
    }
    return 0;
L40:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L45: */
	ewt[i__] = rtol[i__] * (r__1 = ycur[i__], dabs(r__1)) + atol[i__];
    }
    return 0;
/* ----------------------- END OF SUBROUTINE SEWSET ---------------------- */
} /* sewset_ */

