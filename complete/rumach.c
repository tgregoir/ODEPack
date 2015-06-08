/* rumach.f -- translated by f2c (version 20100827).
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

/* Table of constant values */

static real c_b3 = 1.f;

/* DECK RUMACH */
doublereal rumach_(void)
{
    /* System generated locals */
    real ret_val;

    /* Local variables */
    static real u, comp;
    extern /* Subroutine */ int rumsum_(real *, real *, real *);

/* ***BEGIN PROLOGUE  RUMACH */
/* ***PURPOSE  Compute the unit roundoff of the machine. */
/* ***CATEGORY  R1 */
/* ***TYPE      SINGLE PRECISION (RUMACH-S, DUMACH-D) */
/* ***KEYWORDS  MACHINE CONSTANTS */
/* ***AUTHOR  Hindmarsh, Alan C., (LLNL) */
/* ***DESCRIPTION */
/* *Usage: */
/*        REAL  A, RUMACH */
/*        A = RUMACH() */

/* *Function Return Values: */
/*     A : the unit roundoff of the machine. */

/* *Description: */
/*     The unit roundoff is defined as the smallest positive machine */
/*     number u such that  1.0 + u .ne. 1.0.  This is computed by RUMACH */
/*     in a machine-independent manner. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  RUMSUM */
/* ***REVISION HISTORY  (YYYYMMDD) */
/*   19930216  DATE WRITTEN */
/*   19930818  Added SLATEC-format prologue.  (FNF) */
/*   20030707  Added RUMSUM to force normal storage of COMP.  (ACH) */
/* ***END PROLOGUE  RUMACH */

/* ***FIRST EXECUTABLE STATEMENT  RUMACH */
    u = 1.f;
L10:
    u *= .5f;
    rumsum_(&c_b3, &u, &comp);
    if (comp != 1.f) {
	goto L10;
    }
    ret_val = u * 2.f;
    return ret_val;
/* ----------------------- End of Function RUMACH ------------------------ */
} /* rumach_ */

