#ifdef __cplusplus
extern "C" {
#endif

#include "f2c.h"

/* DECK RUMACH */ __device__
doublereal rumach_(void)
{
    /* Table of constant values */

    real c_b3 = 1.f;

    /* System generated locals */
    real ret_val;

    /* Local variables */
    real u, comp;
    extern /* Subroutine */ __device__ int rumsum_(real *, real *, real *);

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

#ifdef __cplusplus
}
#endif
