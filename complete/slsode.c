/* slsode.f -- translated by f2c (version 20100827).
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
    integer init, mxstep, mxhnil, nhnil, nslast, nyh, iowns[6], icf, ierpj, 
	    iersl, jcur, jstart, kflag, l, lyh, lewt, lacor, lsavf, lwm, liwm,
	     meth, miter, maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, 
	    nqu;
} sls001_;

#define sls001_1 sls001_

/* Table of constant values */

static integer c__0 = 0;
static integer c__50 = 50;
static integer c__101 = 101;
static real c_b76 = 0.f;
static integer c__60 = 60;
static integer c__2 = 2;
static integer c__102 = 102;
static integer c__1 = 1;
static integer c__201 = 201;
static integer c__202 = 202;
static integer c__203 = 203;
static integer c__204 = 204;
static integer c__205 = 205;
static integer c__30 = 30;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c__5 = 5;
static integer c__6 = 6;
static integer c__7 = 7;
static integer c__8 = 8;
static integer c__9 = 9;
static integer c__10 = 10;
static integer c__11 = 11;
static integer c__12 = 12;
static integer c__13 = 13;
static integer c__40 = 40;
static integer c__14 = 14;
static integer c__15 = 15;
static integer c__16 = 16;
static integer c__17 = 17;
static integer c__18 = 18;
static integer c__19 = 19;
static integer c__20 = 20;
static integer c__21 = 21;
static integer c__22 = 22;
static integer c__23 = 23;
static integer c__24 = 24;
static integer c__25 = 25;
static integer c__26 = 26;
static integer c__27 = 27;
static integer c__303 = 303;

/* DECK SLSODE */
/* Subroutine */ int slsode_(S_fp f, integer *neq, real *y, real *t, real *
	tout, integer *itol, real *rtol, real *atol, integer *itask, integer *
	istate, integer *iopt, real *rwork, integer *lrw, integer *iwork, 
	integer *liw, U_fp jac, integer *mf)
{
    /* Initialized data */

    static integer mord[2] = { 12,5 };
    static integer mxstp0 = 500;
    static integer mxhnl0 = 10;

    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2;

    /* Builtin functions */
    double sqrt(doublereal), r_sign(real *, real *);

    /* Local variables */
    static integer i__;
    static real h0;
    static integer i1, i2;
    static real w0;
    static integer ml;
    static real rh;
    static integer mu;
    static real tp;
    static integer lf0;
    static real big;
    static integer kgo;
    static real ayi;
    static char msg[80];
    static real hmx, tol, sum, hmax;
    static logical ihit;
    static real ewti, size;
    static integer iflag;
    static real atoli;
    static integer leniw, lenwm, imxer;
    static real tcrit;
    static integer lenrw;
    static real tdist, rtoli, tolsf, tnext;
    extern doublereal rumach_(void);
    extern /* Subroutine */ int sstode_(integer *, real *, real *, integer *, 
	    real *, real *, real *, real *, real *, integer *, S_fp, U_fp, 
	    U_fp, U_fp);
    extern /* Subroutine */ int sprepj_();
    extern /* Subroutine */ int sewset_(integer *, integer *, real *, real *, 
	    real *, real *), sintdy_(real *, integer *, real *, integer *, 
	    real *, integer *);
    extern doublereal svnorm_(integer *, real *, real *);
    extern /* Subroutine */ int ssolsy_();
    extern /* Subroutine */ int xerrwv_(char *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, real *, real *, 
	    ftnlen);

/* ***BEGIN PROLOGUE  SLSODE */
/* ***PURPOSE  Livermore Solver for Ordinary Differential Equations. */
/*            SLSODE solves the initial-value problem for stiff or */
/*            nonstiff systems of first-order ODE's, */
/*               dy/dt = f(t,y),   or, in component form, */
/*               dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(N)),  i=1,...,N. */
/* ***CATEGORY  I1A */
/* ***TYPE      SINGLE PRECISION (SLSODE-S, DLSODE-D) */
/* ***KEYWORDS  ORDINARY DIFFERENTIAL EQUATIONS, INITIAL VALUE PROBLEM, */
/*             STIFF, NONSTIFF */
/* ***AUTHOR  Hindmarsh, Alan C., (LLNL) */
/*             Center for Applied Scientific Computing, L-561 */
/*             Lawrence Livermore National Laboratory */
/*             Livermore, CA 94551. */
/* ***DESCRIPTION */

/*     NOTE: The "Usage" and "Arguments" sections treat only a subset of */
/*           available options, in condensed fashion.  The options */
/*           covered and the information supplied will support most */
/*           standard uses of SLSODE. */

/*           For more sophisticated uses, full details on all options are */
/*           given in the concluding section, headed "Long Description." */
/*           A synopsis of the SLSODE Long Description is provided at the */
/*           beginning of that section; general topics covered are: */
/*           - Elements of the call sequence; optional input and output */
/*           - Optional supplemental routines in the SLSODE package */
/*           - internal COMMON block */

/* *Usage: */
/*     Communication between the user and the SLSODE package, for normal */
/*     situations, is summarized here.  This summary describes a subset */
/*     of the available options.  See "Long Description" for complete */
/*     details, including optional communication, nonstandard options, */
/*     and instructions for special situations. */

/*     A sample program is given in the "Examples" section. */

/*     Refer to the argument descriptions for the definitions of the */
/*     quantities that appear in the following sample declarations. */

/*     For MF = 10, */
/*        PARAMETER  (LRW = 20 + 16*NEQ,           LIW = 20) */
/*     For MF = 21 or 22, */
/*        PARAMETER  (LRW = 22 +  9*NEQ + NEQ**2,  LIW = 20 + NEQ) */
/*     For MF = 24 or 25, */
/*        PARAMETER  (LRW = 22 + 10*NEQ + (2*ML+MU)*NEQ, */
/*       *                                         LIW = 20 + NEQ) */

/*        EXTERNAL F, JAC */
/*        INTEGER  NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK(LIW), */
/*       *         LIW, MF */
/*        REAL Y(NEQ), T, TOUT, RTOL, ATOL(ntol), RWORK(LRW) */

/*        CALL SLSODE (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, */
/*       *            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF) */

/* *Arguments: */
/*     F     :EXT    Name of subroutine for right-hand-side vector f. */
/*                   This name must be declared EXTERNAL in calling */
/*                   program.  The form of F must be: */

/*                   SUBROUTINE  F (NEQ, T, Y, YDOT) */
/*                   INTEGER  NEQ */
/*                   REAL T, Y(*), YDOT(*) */

/*                   The inputs are NEQ, T, Y.  F is to set */

/*                   YDOT(i) = f(i,T,Y(1),Y(2),...,Y(NEQ)), */
/*                                                     i = 1, ..., NEQ . */

/*     NEQ   :IN     Number of first-order ODE's. */

/*     Y     :INOUT  Array of values of the y(t) vector, of length NEQ. */
/*                   Input:  For the first call, Y should contain the */
/*                           values of y(t) at t = T. (Y is an input */
/*                           variable only if ISTATE = 1.) */
/*                   Output: On return, Y will contain the values at the */
/*                           new t-value. */

/*     T     :INOUT  Value of the independent variable.  On return it */
/*                   will be the current value of t (normally TOUT). */

/*     TOUT  :IN     Next point where output is desired (.NE. T). */

/*     ITOL  :IN     1 or 2 according as ATOL (below) is a scalar or */
/*                   an array. */

/*     RTOL  :IN     Relative tolerance parameter (scalar). */

/*     ATOL  :IN     Absolute tolerance parameter (scalar or array). */
/*                   If ITOL = 1, ATOL need not be dimensioned. */
/*                   If ITOL = 2, ATOL must be dimensioned at least NEQ. */

/*                   The estimated local error in Y(i) will be controlled */
/*                   so as to be roughly less (in magnitude) than */

/*                   EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or */
/*                   EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2. */

/*                   Thus the local error test passes if, in each */
/*                   component, either the absolute error is less than */
/*                   ATOL (or ATOL(i)), or the relative error is less */
/*                   than RTOL. */

/*                   Use RTOL = 0.0 for pure absolute error control, and */
/*                   use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative */
/*                   error control.  Caution:  Actual (global) errors may */
/*                   exceed these local tolerances, so choose them */
/*                   conservatively. */

/*     ITASK :IN     Flag indicating the task SLSODE is to perform. */
/*                   Use ITASK = 1 for normal computation of output */
/*                   values of y at t = TOUT. */

/*     ISTATE:INOUT  Index used for input and output to specify the state */
/*                   of the calculation. */
/*                   Input: */
/*                    1   This is the first call for a problem. */
/*                    2   This is a subsequent call. */
/*                   Output: */
/*                    1   Nothing was done, as TOUT was equal to T. */
/*                    2   SLSODE was successful (otherwise, negative). */
/*                        Note that ISTATE need not be modified after a */
/*                        successful return. */
/*                   -1   Excess work done on this call (perhaps wrong */
/*                        MF). */
/*                   -2   Excess accuracy requested (tolerances too */
/*                        small). */
/*                   -3   Illegal input detected (see printed message). */
/*                   -4   Repeated error test failures (check all */
/*                        inputs). */
/*                   -5   Repeated convergence failures (perhaps bad */
/*                        Jacobian supplied or wrong choice of MF or */
/*                        tolerances). */
/*                   -6   Error weight became zero during problem */
/*                        (solution component i vanished, and ATOL or */
/*                        ATOL(i) = 0.). */

/*     IOPT  :IN     Flag indicating whether optional inputs are used: */
/*                   0   No. */
/*                   1   Yes.  (See "Optional inputs" under "Long */
/*                       Description," Part 1.) */

/*     RWORK :WORK   Real work array of length at least: */
/*                   20 + 16*NEQ                    for MF = 10, */
/*                   22 +  9*NEQ + NEQ**2           for MF = 21 or 22, */
/*                   22 + 10*NEQ + (2*ML + MU)*NEQ  for MF = 24 or 25. */

/*     LRW   :IN     Declared length of RWORK (in user's DIMENSION */
/*                   statement). */

/*     IWORK :WORK   Integer work array of length at least: */
/*                   20        for MF = 10, */
/*                   20 + NEQ  for MF = 21, 22, 24, or 25. */

/*                   If MF = 24 or 25, input in IWORK(1),IWORK(2) the */
/*                   lower and upper Jacobian half-bandwidths ML,MU. */

/*                   On return, IWORK contains information that may be */
/*                   of interest to the user: */

/*            Name   Location   Meaning */
/*            -----  ---------  ----------------------------------------- */
/*            NST    IWORK(11)  Number of steps taken for the problem so */
/*                              far. */
/*            NFE    IWORK(12)  Number of f evaluations for the problem */
/*                              so far. */
/*            NJE    IWORK(13)  Number of Jacobian evaluations (and of */
/*                              matrix LU decompositions) for the problem */
/*                              so far. */
/*            NQU    IWORK(14)  Method order last used (successfully). */
/*            LENRW  IWORK(17)  Length of RWORK actually required.  This */
/*                              is defined on normal returns and on an */
/*                              illegal input return for insufficient */
/*                              storage. */
/*            LENIW  IWORK(18)  Length of IWORK actually required.  This */
/*                              is defined on normal returns and on an */
/*                              illegal input return for insufficient */
/*                              storage. */

/*     LIW   :IN     Declared length of IWORK (in user's DIMENSION */
/*                   statement). */

/*     JAC   :EXT    Name of subroutine for Jacobian matrix (MF = */
/*                   21 or 24).  If used, this name must be declared */
/*                   EXTERNAL in calling program.  If not used, pass a */
/*                   dummy name.  The form of JAC must be: */

/*                   SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD) */
/*                   INTEGER NEQ, ML, MU, NROWPD */
/*                   REAL T, Y(*), PD(NROWPD,*) */

/*                   See item c, under "Description" below for more */
/*                   information about JAC. */

/*     MF    :IN     Method flag.  Standard values are: */
/*                   10  Nonstiff (Adams) method, no Jacobian used. */
/*                   21  Stiff (BDF) method, user-supplied full Jacobian. */
/*                   22  Stiff method, internally generated full */
/*                       Jacobian. */
/*                   24  Stiff method, user-supplied banded Jacobian. */
/*                   25  Stiff method, internally generated banded */
/*                       Jacobian. */

/* *Description: */
/*     SLSODE solves the initial value problem for stiff or nonstiff */
/*     systems of first-order ODE's, */

/*        dy/dt = f(t,y) , */

/*     or, in component form, */

/*        dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ)) */
/*                                                  (i = 1, ..., NEQ) . */

/*     SLSODE is a package based on the GEAR and GEARB packages, and on */
/*     the October 23, 1978, version of the tentative ODEPACK user */
/*     interface standard, with minor modifications. */

/*     The steps in solving such a problem are as follows. */

/*     a. First write a subroutine of the form */

/*           SUBROUTINE  F (NEQ, T, Y, YDOT) */
/*           INTEGER  NEQ */
/*           REAL T, Y(*), YDOT(*) */

/*        which supplies the vector function f by loading YDOT(i) with */
/*        f(i). */

/*     b. Next determine (or guess) whether or not the problem is stiff. */
/*        Stiffness occurs when the Jacobian matrix df/dy has an */
/*        eigenvalue whose real part is negative and large in magnitude */
/*        compared to the reciprocal of the t span of interest.  If the */
/*        problem is nonstiff, use method flag MF = 10.  If it is stiff, */
/*        there are four standard choices for MF, and SLSODE requires the */
/*        Jacobian matrix in some form.  This matrix is regarded either */
/*        as full (MF = 21 or 22), or banded (MF = 24 or 25).  In the */
/*        banded case, SLSODE requires two half-bandwidth parameters ML */
/*        and MU. These are, respectively, the widths of the lower and */
/*        upper parts of the band, excluding the main diagonal.  Thus the */
/*        band consists of the locations (i,j) with */

/*           i - ML <= j <= i + MU , */

/*        and the full bandwidth is ML + MU + 1 . */

/*     c. If the problem is stiff, you are encouraged to supply the */
/*        Jacobian directly (MF = 21 or 24), but if this is not feasible, */
/*        SLSODE will compute it internally by difference quotients (MF = */
/*        22 or 25).  If you are supplying the Jacobian, write a */
/*        subroutine of the form */

/*           SUBROUTINE  JAC (NEQ, T, Y, ML, MU, PD, NROWPD) */
/*           INTEGER  NEQ, ML, MU, NRWOPD */
/*           REAL T, Y(*), PD(NROWPD,*) */

/*        which provides df/dy by loading PD as follows: */
/*        - For a full Jacobian (MF = 21), load PD(i,j) with df(i)/dy(j), */
/*          the partial derivative of f(i) with respect to y(j).  (Ignore */
/*          the ML and MU arguments in this case.) */
/*        - For a banded Jacobian (MF = 24), load PD(i-j+MU+1,j) with */
/*          df(i)/dy(j); i.e., load the diagonal lines of df/dy into the */
/*          rows of PD from the top down. */
/*        - In either case, only nonzero elements need be loaded. */

/*     d. Write a main program that calls subroutine SLSODE once for each */
/*        point at which answers are desired.  This should also provide */
/*        for possible use of logical unit 6 for output of error messages */
/*        by SLSODE. */

/*        Before the first call to SLSODE, set ISTATE = 1, set Y and T to */
/*        the initial values, and set TOUT to the first output point.  To */
/*        continue the integration after a successful return, simply */
/*        reset TOUT and call SLSODE again.  No other parameters need be */
/*        reset. */

/* *Examples: */
/*     The following is a simple example problem, with the coding needed */
/*     for its solution by SLSODE. The problem is from chemical kinetics, */
/*     and consists of the following three rate equations: */

/*        dy1/dt = -.04*y1 + 1.E4*y2*y3 */
/*        dy2/dt = .04*y1 - 1.E4*y2*y3 - 3.E7*y2**2 */
/*        dy3/dt = 3.E7*y2**2 */

/*     on the interval from t = 0.0 to t = 4.E10, with initial conditions */
/*     y1 = 1.0, y2 = y3 = 0. The problem is stiff. */

/*     The following coding solves this problem with SLSODE, using */
/*     MF = 21 and printing results at t = .4, 4., ..., 4.E10.  It uses */
/*     ITOL = 2 and ATOL much smaller for y2 than for y1 or y3 because y2 */
/*     has much smaller values.  At the end of the run, statistical */
/*     quantities of interest are printed. */

/*        EXTERNAL  FEX, JEX */
/*        INTEGER  IOPT, IOUT, ISTATE, ITASK, ITOL, IWORK(23), LIW, LRW, */
/*       *         MF, NEQ */
/*        REAL  ATOL(3), RTOL, RWORK(58), T, TOUT, Y(3) */
/*        NEQ = 3 */
/*        Y(1) = 1. */
/*        Y(2) = 0. */
/*        Y(3) = 0. */
/*        T = 0. */
/*        TOUT = .4 */
/*        ITOL = 2 */
/*        RTOL = 1.E-4 */
/*        ATOL(1) = 1.E-6 */
/*        ATOL(2) = 1.E-10 */
/*        ATOL(3) = 1.E-6 */
/*        ITASK = 1 */
/*        ISTATE = 1 */
/*        IOPT = 0 */
/*        LRW = 58 */
/*        LIW = 23 */
/*        MF = 21 */
/*        DO 40 IOUT = 1,12 */
/*          CALL SLSODE (FEX, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, */
/*       *               ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JEX, MF) */
/*          WRITE(6,20)  T, Y(1), Y(2), Y(3) */
/*    20    FORMAT(' At t =',E12.4,'   y =',3E14.6) */
/*          IF (ISTATE .LT. 0)  GO TO 80 */
/*    40    TOUT = TOUT*10. */
/*        WRITE(6,60)  IWORK(11), IWORK(12), IWORK(13) */
/*    60  FORMAT(/' No. steps =',i4,',  No. f-s =',i4,',  No. J-s =',i4) */
/*        STOP */
/*    80  WRITE(6,90)  ISTATE */
/*    90  FORMAT(///' Error halt.. ISTATE =',I3) */
/*        STOP */
/*        END */

/*        SUBROUTINE  FEX (NEQ, T, Y, YDOT) */
/*        INTEGER  NEQ */
/*        REAL  T, Y(3), YDOT(3) */
/*        YDOT(1) = -.04*Y(1) + 1.E4*Y(2)*Y(3) */
/*        YDOT(3) = 3.E7*Y(2)*Y(2) */
/*        YDOT(2) = -YDOT(1) - YDOT(3) */
/*        RETURN */
/*        END */

/*        SUBROUTINE  JEX (NEQ, T, Y, ML, MU, PD, NRPD) */
/*        INTEGER  NEQ, ML, MU, NRPD */
/*        REAL  T, Y(3), PD(NRPD,3) */
/*        PD(1,1) = -.04 */
/*        PD(1,2) = 1.E4*Y(3) */
/*        PD(1,3) = 1.E4*Y(2) */
/*        PD(2,1) = .04 */
/*        PD(2,3) = -PD(1,3) */
/*        PD(3,2) = 6.E7*Y(2) */
/*        PD(2,2) = -PD(1,2) - PD(3,2) */
/*        RETURN */
/*        END */

/*     The output from this program (on a Cray-1 in single precision) */
/*     is as follows. */

/*     At t =  4.0000e-01   y =  9.851726e-01  3.386406e-05  1.479357e-02 */
/*     At t =  4.0000e+00   y =  9.055142e-01  2.240418e-05  9.446344e-02 */
/*     At t =  4.0000e+01   y =  7.158050e-01  9.184616e-06  2.841858e-01 */
/*     At t =  4.0000e+02   y =  4.504846e-01  3.222434e-06  5.495122e-01 */
/*     At t =  4.0000e+03   y =  1.831701e-01  8.940379e-07  8.168290e-01 */
/*     At t =  4.0000e+04   y =  3.897016e-02  1.621193e-07  9.610297e-01 */
/*     At t =  4.0000e+05   y =  4.935213e-03  1.983756e-08  9.950648e-01 */
/*     At t =  4.0000e+06   y =  5.159269e-04  2.064759e-09  9.994841e-01 */
/*     At t =  4.0000e+07   y =  5.306413e-05  2.122677e-10  9.999469e-01 */
/*     At t =  4.0000e+08   y =  5.494530e-06  2.197825e-11  9.999945e-01 */
/*     At t =  4.0000e+09   y =  5.129458e-07  2.051784e-12  9.999995e-01 */
/*     At t =  4.0000e+10   y = -7.170603e-08 -2.868241e-13  1.000000e+00 */

/*     No. steps = 330,  No. f-s = 405,  No. J-s = 69 */

/* *Accuracy: */
/*     The accuracy of the solution depends on the choice of tolerances */
/*     RTOL and ATOL.  Actual (global) errors may exceed these local */
/*     tolerances, so choose them conservatively. */

/* *Cautions: */
/*     The work arrays should not be altered between calls to SLSODE for */
/*     the same problem, except possibly for the conditional and optional */
/*     inputs. */

/* *Portability: */
/*     Since NEQ is dimensioned inside SLSODE, some compilers may object */
/*     to a call to SLSODE with NEQ a scalar variable.  In this event, */
/*     use DIMENSION NEQ(1).  Similar remarks apply to RTOL and ATOL. */

/*     Note to Cray users: */
/*     For maximum efficiency, use the CFT77 compiler.  Appropriate */
/*     compiler optimization directives have been inserted for CFT77. */

/* *Reference: */
/*     Alan C. Hindmarsh, "ODEPACK, A Systematized Collection of ODE */
/*     Solvers," in Scientific Computing, R. S. Stepleman, et al., Eds. */
/*     (North-Holland, Amsterdam, 1983), pp. 55-64. */

/* *Long Description: */
/*     The following complete description of the user interface to */
/*     SLSODE consists of four parts: */

/*     1.  The call sequence to subroutine SLSODE, which is a driver */
/*         routine for the solver.  This includes descriptions of both */
/*         the call sequence arguments and user-supplied routines. */
/*         Following these descriptions is a description of optional */
/*         inputs available through the call sequence, and then a */
/*         description of optional outputs in the work arrays. */

/*     2.  Descriptions of other routines in the SLSODE package that may */
/*         be (optionally) called by the user.  These provide the ability */
/*         to alter error message handling, save and restore the internal */
/*         COMMON, and obtain specified derivatives of the solution y(t). */

/*     3.  Descriptions of COMMON block to be declared in overlay or */
/*         similar environments, or to be saved when doing an interrupt */
/*         of the problem and continued solution later. */

/*     4.  Description of two routines in the SLSODE package, either of */
/*         which the user may replace with his own version, if desired. */
/*         These relate to the measurement of errors. */


/*                         Part 1.  Call Sequence */
/*                         ---------------------- */

/*     Arguments */
/*     --------- */
/*     The call sequence parameters used for input only are */

/*        F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW, JAC, MF, */

/*     and those used for both input and output are */

/*        Y, T, ISTATE. */

/*     The work arrays RWORK and IWORK are also used for conditional and */
/*     optional inputs and optional outputs.  (The term output here */
/*     refers to the return from subroutine SLSODE to the user's calling */
/*     program.) */

/*     The legality of input parameters will be thoroughly checked on the */
/*     initial call for the problem, but not checked thereafter unless a */
/*     change in input parameters is flagged by ISTATE = 3 on input. */

/*     The descriptions of the call arguments are as follows. */

/*     F        The name of the user-supplied subroutine defining the ODE */
/*              system.  The system must be put in the first-order form */
/*              dy/dt = f(t,y), where f is a vector-valued function of */
/*              the scalar t and the vector y. Subroutine F is to compute */
/*              the function f. It is to have the form */

/*                 SUBROUTINE F (NEQ, T, Y, YDOT) */
/*                 REAL T, Y(*), YDOT(*) */

/*              where NEQ, T, and Y are input, and the array YDOT = */
/*              f(T,Y) is output.  Y and YDOT are arrays of length NEQ. */
/*              Subroutine F should not alter Y(1),...,Y(NEQ).  F must be */
/*              declared EXTERNAL in the calling program. */

/*              Subroutine F may access user-defined quantities in */
/*              NEQ(2),... and/or in Y(NEQ(1)+1),..., if NEQ is an array */
/*              (dimensioned in F) and/or Y has length exceeding NEQ(1). */
/*              See the descriptions of NEQ and Y below. */

/*              If quantities computed in the F routine are needed */
/*              externally to SLSODE, an extra call to F should be made */
/*              for this purpose, for consistent and accurate results. */
/*              If only the derivative dy/dt is needed, use SINTDY */
/*              instead. */

/*     NEQ      The size of the ODE system (number of first-order */
/*              ordinary differential equations).  Used only for input. */
/*              NEQ may be decreased, but not increased, during the */
/*              problem.  If NEQ is decreased (with ISTATE = 3 on input), */
/*              the remaining components of Y should be left undisturbed, */
/*              if these are to be accessed in F and/or JAC. */

/*              Normally, NEQ is a scalar, and it is generally referred */
/*              to as a scalar in this user interface description. */
/*              However, NEQ may be an array, with NEQ(1) set to the */
/*              system size.  (The SLSODE package accesses only NEQ(1).) */
/*              In either case, this parameter is passed as the NEQ */
/*              argument in all calls to F and JAC.  Hence, if it is an */
/*              array, locations NEQ(2),... may be used to store other */
/*              integer data and pass it to F and/or JAC.  Subroutines */
/*              F and/or JAC must include NEQ in a DIMENSION statement */
/*              in that case. */

/*     Y        A real array for the vector of dependent variables, of */
/*              length NEQ or more.  Used for both input and output on */
/*              the first call (ISTATE = 1), and only for output on */
/*              other calls.  On the first call, Y must contain the */
/*              vector of initial values.  On output, Y contains the */
/*              computed solution vector, evaluated at T. If desired, */
/*              the Y array may be used for other purposes between */
/*              calls to the solver. */

/*              This array is passed as the Y argument in all calls to F */
/*              and JAC.  Hence its length may exceed NEQ, and locations */
/*              Y(NEQ+1),... may be used to store other real data and */
/*              pass it to F and/or JAC.  (The SLSODE package accesses */
/*              only Y(1),...,Y(NEQ).) */

/*     T        The independent variable.  On input, T is used only on */
/*              the first call, as the initial point of the integration. */
/*              On output, after each call, T is the value at which a */
/*              computed solution Y is evaluated (usually the same as */
/*              TOUT).  On an error return, T is the farthest point */
/*              reached. */

/*     TOUT     The next value of T at which a computed solution is */
/*              desired.  Used only for input. */

/*              When starting the problem (ISTATE = 1), TOUT may be equal */
/*              to T for one call, then should not equal T for the next */
/*              call.  For the initial T, an input value of TOUT .NE. T */
/*              is used in order to determine the direction of the */
/*              integration (i.e., the algebraic sign of the step sizes) */
/*              and the rough scale of the problem.  Integration in */
/*              either direction (forward or backward in T) is permitted. */

/*              If ITASK = 2 or 5 (one-step modes), TOUT is ignored */
/*              after the first call (i.e., the first call with */
/*              TOUT .NE. T).  Otherwise, TOUT is required on every call. */

/*              If ITASK = 1, 3, or 4, the values of TOUT need not be */
/*              monotone, but a value of TOUT which backs up is limited */
/*              to the current internal T interval, whose endpoints are */
/*              TCUR - HU and TCUR.  (See "Optional Outputs" below for */
/*              TCUR and HU.) */


/*     ITOL     An indicator for the type of error control.  See */
/*              description below under ATOL.  Used only for input. */

/*     RTOL     A relative error tolerance parameter, either a scalar or */
/*              an array of length NEQ.  See description below under */
/*              ATOL.  Input only. */

/*     ATOL     An absolute error tolerance parameter, either a scalar or */
/*              an array of length NEQ.  Input only. */

/*              The input parameters ITOL, RTOL, and ATOL determine the */
/*              error control performed by the solver.  The solver will */
/*              control the vector e = (e(i)) of estimated local errors */
/*              in Y, according to an inequality of the form */

/*                 rms-norm of ( e(i)/EWT(i) ) <= 1, */

/*              where */

/*                 EWT(i) = RTOL(i)*ABS(Y(i)) + ATOL(i), */

/*              and the rms-norm (root-mean-square norm) here is */

/*                 rms-norm(v) = SQRT(sum v(i)**2 / NEQ). */

/*              Here EWT = (EWT(i)) is a vector of weights which must */
/*              always be positive, and the values of RTOL and ATOL */
/*              should all be nonnegative.  The following table gives the */
/*              types (scalar/array) of RTOL and ATOL, and the */
/*              corresponding form of EWT(i). */

/*              ITOL    RTOL      ATOL      EWT(i) */
/*              ----    ------    ------    ----------------------------- */
/*              1       scalar    scalar    RTOL*ABS(Y(i)) + ATOL */
/*              2       scalar    array     RTOL*ABS(Y(i)) + ATOL(i) */
/*              3       array     scalar    RTOL(i)*ABS(Y(i)) + ATOL */
/*              4       array     array     RTOL(i)*ABS(Y(i)) + ATOL(i) */

/*              When either of these parameters is a scalar, it need not */
/*              be dimensioned in the user's calling program. */

/*              If none of the above choices (with ITOL, RTOL, and ATOL */
/*              fixed throughout the problem) is suitable, more general */
/*              error controls can be obtained by substituting */
/*              user-supplied routines for the setting of EWT and/or for */
/*              the norm calculation.  See Part 4 below. */

/*              If global errors are to be estimated by making a repeated */
/*              run on the same problem with smaller tolerances, then all */
/*              components of RTOL and ATOL (i.e., of EWT) should be */
/*              scaled down uniformly. */

/*     ITASK    An index specifying the task to be performed.  Input */
/*              only.  ITASK has the following values and meanings: */
/*              1   Normal computation of output values of y(t) at */
/*                  t = TOUT (by overshooting and interpolating). */
/*              2   Take one step only and return. */
/*              3   Stop at the first internal mesh point at or beyond */
/*                  t = TOUT and return. */
/*              4   Normal computation of output values of y(t) at */
/*                  t = TOUT but without overshooting t = TCRIT.  TCRIT */
/*                  must be input as RWORK(1).  TCRIT may be equal to or */
/*                  beyond TOUT, but not behind it in the direction of */
/*                  integration.  This option is useful if the problem */
/*                  has a singularity at or beyond t = TCRIT. */
/*              5   Take one step, without passing TCRIT, and return. */
/*                  TCRIT must be input as RWORK(1). */

/*              Note:  If ITASK = 4 or 5 and the solver reaches TCRIT */
/*              (within roundoff), it will return T = TCRIT (exactly) to */
/*              indicate this (unless ITASK = 4 and TOUT comes before */
/*              TCRIT, in which case answers at T = TOUT are returned */
/*              first). */

/*     ISTATE   An index used for input and output to specify the state */
/*              of the calculation. */

/*              On input, the values of ISTATE are as follows: */
/*              1   This is the first call for the problem */
/*                  (initializations will be done).  See "Note" below. */
/*              2   This is not the first call, and the calculation is to */
/*                  continue normally, with no change in any input */
/*                  parameters except possibly TOUT and ITASK.  (If ITOL, */
/*                  RTOL, and/or ATOL are changed between calls with */
/*                  ISTATE = 2, the new values will be used but not */
/*                  tested for legality.) */
/*              3   This is not the first call, and the calculation is to */
/*                  continue normally, but with a change in input */
/*                  parameters other than TOUT and ITASK.  Changes are */
/*                  allowed in NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, MF, */
/*                  ML, MU, and any of the optional inputs except H0. */
/*                  (See IWORK description for ML and MU.) */

/*              Note:  A preliminary call with TOUT = T is not counted as */
/*              a first call here, as no initialization or checking of */
/*              input is done.  (Such a call is sometimes useful for the */
/*              purpose of outputting the initial conditions.)  Thus the */
/*              first call for which TOUT .NE. T requires ISTATE = 1 on */
/*              input. */

/*              On output, ISTATE has the following values and meanings: */
/*               1  Nothing was done, as TOUT was equal to T with */
/*                  ISTATE = 1 on input. */
/*               2  The integration was performed successfully. */
/*              -1  An excessive amount of work (more than MXSTEP steps) */
/*                  was done on this call, before completing the */
/*                  requested task, but the integration was otherwise */
/*                  successful as far as T. (MXSTEP is an optional input */
/*                  and is normally 500.)  To continue, the user may */
/*                  simply reset ISTATE to a value >1 and call again (the */
/*                  excess work step counter will be reset to 0).  In */
/*                  addition, the user may increase MXSTEP to avoid this */
/*                  error return; see "Optional Inputs" below. */
/*              -2  Too much accuracy was requested for the precision of */
/*                  the machine being used.  This was detected before */
/*                  completing the requested task, but the integration */
/*                  was successful as far as T. To continue, the */
/*                  tolerance parameters must be reset, and ISTATE must */
/*                  be set to 3. The optional output TOLSF may be used */
/*                  for this purpose.  (Note:  If this condition is */
/*                  detected before taking any steps, then an illegal */
/*                  input return (ISTATE = -3) occurs instead.) */
/*              -3  Illegal input was detected, before taking any */
/*                  integration steps.  See written message for details. */
/*                  (Note:  If the solver detects an infinite loop of */
/*                  calls to the solver with illegal input, it will cause */
/*                  the run to stop.) */
/*              -4  There were repeated error-test failures on one */
/*                  attempted step, before completing the requested task, */
/*                  but the integration was successful as far as T.  The */
/*                  problem may have a singularity, or the input may be */
/*                  inappropriate. */
/*              -5  There were repeated convergence-test failures on one */
/*                  attempted step, before completing the requested task, */
/*                  but the integration was successful as far as T. This */
/*                  may be caused by an inaccurate Jacobian matrix, if */
/*                  one is being used. */
/*              -6  EWT(i) became zero for some i during the integration. */
/*                  Pure relative error control (ATOL(i)=0.0) was */
/*                  requested on a variable which has now vanished.  The */
/*                  integration was successful as far as T. */

/*              Note:  Since the normal output value of ISTATE is 2, it */
/*              does not need to be reset for normal continuation.  Also, */
/*              since a negative input value of ISTATE will be regarded */
/*              as illegal, a negative output value requires the user to */
/*              change it, and possibly other inputs, before calling the */
/*              solver again. */

/*     IOPT     An integer flag to specify whether any optional inputs */
/*              are being used on this call.  Input only.  The optional */
/*              inputs are listed under a separate heading below. */
/*              0   No optional inputs are being used.  Default values */
/*                  will be used in all cases. */
/*              1   One or more optional inputs are being used. */

/*     RWORK    A real working array (single precision).  The length of */
/*              RWORK must be at least */

/*                 20 + NYH*(MAXORD + 1) + 3*NEQ + LWM */

/*              where */
/*                 NYH = the initial value of NEQ, */
/*              MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a */
/*                       smaller value is given as an optional input), */
/*                 LWM = 0           if MITER = 0, */
/*                 LWM = NEQ**2 + 2  if MITER = 1 or 2, */
/*                 LWM = NEQ + 2     if MITER = 3, and */
/*                 LWM = (2*ML + MU + 1)*NEQ + 2 */
/*                                   if MITER = 4 or 5. */
/*              (See the MF description below for METH and MITER.) */

/*              Thus if MAXORD has its default value and NEQ is constant, */
/*              this length is: */
/*              20 + 16*NEQ                    for MF = 10, */
/*              22 + 16*NEQ + NEQ**2           for MF = 11 or 12, */
/*              22 + 17*NEQ                    for MF = 13, */
/*              22 + 17*NEQ + (2*ML + MU)*NEQ  for MF = 14 or 15, */
/*              20 +  9*NEQ                    for MF = 20, */
/*              22 +  9*NEQ + NEQ**2           for MF = 21 or 22, */
/*              22 + 10*NEQ                    for MF = 23, */
/*              22 + 10*NEQ + (2*ML + MU)*NEQ  for MF = 24 or 25. */

/*              The first 20 words of RWORK are reserved for conditional */
/*              and optional inputs and optional outputs. */

/*              The following word in RWORK is a conditional input: */
/*              RWORK(1) = TCRIT, the critical value of t which the */
/*                         solver is not to overshoot.  Required if ITASK */
/*                         is 4 or 5, and ignored otherwise.  See ITASK. */

/*     LRW      The length of the array RWORK, as declared by the user. */
/*              (This will be checked by the solver.) */

/*     IWORK    An integer work array.  Its length must be at least */
/*              20       if MITER = 0 or 3 (MF = 10, 13, 20, 23), or */
/*              20 + NEQ otherwise (MF = 11, 12, 14, 15, 21, 22, 24, 25). */
/*              (See the MF description below for MITER.)  The first few */
/*              words of IWORK are used for conditional and optional */
/*              inputs and optional outputs. */

/*              The following two words in IWORK are conditional inputs: */
/*              IWORK(1) = ML   These are the lower and upper half- */
/*              IWORK(2) = MU   bandwidths, respectively, of the banded */
/*                              Jacobian, excluding the main diagonal. */
/*                         The band is defined by the matrix locations */
/*                         (i,j) with i - ML <= j <= i + MU. ML and MU */
/*                         must satisfy 0 <= ML,MU <= NEQ - 1. These are */
/*                         required if MITER is 4 or 5, and ignored */
/*                         otherwise.  ML and MU may in fact be the band */
/*                         parameters for a matrix to which df/dy is only */
/*                         approximately equal. */

/*     LIW      The length of the array IWORK, as declared by the user. */
/*              (This will be checked by the solver.) */

/*     Note:  The work arrays must not be altered between calls to SLSODE */
/*     for the same problem, except possibly for the conditional and */
/*     optional inputs, and except for the last 3*NEQ words of RWORK. */
/*     The latter space is used for internal scratch space, and so is */
/*     available for use by the user outside SLSODE between calls, if */
/*     desired (but not for use by F or JAC). */

/*     JAC      The name of the user-supplied routine (MITER = 1 or 4) to */
/*              compute the Jacobian matrix, df/dy, as a function of the */
/*              scalar t and the vector y.  (See the MF description below */
/*              for MITER.)  It is to have the form */

/*                 SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD) */
/*                 REAL T, Y(*), PD(NROWPD,*) */

/*              where NEQ, T, Y, ML, MU, and NROWPD are input and the */
/*              array PD is to be loaded with partial derivatives */
/*              (elements of the Jacobian matrix) on output.  PD must be */
/*              given a first dimension of NROWPD.  T and Y have the same */
/*              meaning as in subroutine F. */

/*              In the full matrix case (MITER = 1), ML and MU are */
/*              ignored, and the Jacobian is to be loaded into PD in */
/*              columnwise manner, with df(i)/dy(j) loaded into PD(i,j). */

/*              In the band matrix case (MITER = 4), the elements within */
/*              the band are to be loaded into PD in columnwise manner, */
/*              with diagonal lines of df/dy loaded into the rows of PD. */
/*              Thus df(i)/dy(j) is to be loaded into PD(i-j+MU+1,j).  ML */
/*              and MU are the half-bandwidth parameters (see IWORK). */
/*              The locations in PD in the two triangular areas which */
/*              correspond to nonexistent matrix elements can be ignored */
/*              or loaded arbitrarily, as they are overwritten by SLSODE. */

/*              JAC need not provide df/dy exactly. A crude approximation */
/*              (possibly with a smaller bandwidth) will do. */

/*              In either case, PD is preset to zero by the solver, so */
/*              that only the nonzero elements need be loaded by JAC. */
/*              Each call to JAC is preceded by a call to F with the same */
/*              arguments NEQ, T, and Y. Thus to gain some efficiency, */
/*              intermediate quantities shared by both calculations may */
/*              be saved in a user COMMON block by F and not recomputed */
/*              by JAC, if desired.  Also, JAC may alter the Y array, if */
/*              desired.  JAC must be declared EXTERNAL in the calling */
/*              program. */

/*              Subroutine JAC may access user-defined quantities in */
/*              NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array */
/*              (dimensioned in JAC) and/or Y has length exceeding */
/*              NEQ(1).  See the descriptions of NEQ and Y above. */

/*     MF       The method flag.  Used only for input.  The legal values */
/*              of MF are 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24, */
/*              and 25.  MF has decimal digits METH and MITER: */
/*                 MF = 10*METH + MITER . */

/*              METH indicates the basic linear multistep method: */
/*              1   Implicit Adams method. */
/*              2   Method based on backward differentiation formulas */
/*                  (BDF's). */

/*              MITER indicates the corrector iteration method: */
/*              0   Functional iteration (no Jacobian matrix is */
/*                  involved). */
/*              1   Chord iteration with a user-supplied full (NEQ by */
/*                  NEQ) Jacobian. */
/*              2   Chord iteration with an internally generated */
/*                  (difference quotient) full Jacobian (using NEQ */
/*                  extra calls to F per df/dy value). */
/*              3   Chord iteration with an internally generated */
/*                  diagonal Jacobian approximation (using one extra call */
/*                  to F per df/dy evaluation). */
/*              4   Chord iteration with a user-supplied banded Jacobian. */
/*              5   Chord iteration with an internally generated banded */
/*                  Jacobian (using ML + MU + 1 extra calls to F per */
/*                  df/dy evaluation). */

/*              If MITER = 1 or 4, the user must supply a subroutine JAC */
/*              (the name is arbitrary) as described above under JAC. */
/*              For other values of MITER, a dummy argument can be used. */

/*     Optional Inputs */
/*     --------------- */
/*     The following is a list of the optional inputs provided for in the */
/*     call sequence.  (See also Part 2.)  For each such input variable, */
/*     this table lists its name as used in this documentation, its */
/*     location in the call sequence, its meaning, and the default value. */
/*     The use of any of these inputs requires IOPT = 1, and in that case */
/*     all of these inputs are examined.  A value of zero for any of */
/*     these optional inputs will cause the default value to be used. */
/*     Thus to use a subset of the optional inputs, simply preload */
/*     locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, */
/*     and then set those of interest to nonzero values. */

/*     Name    Location   Meaning and default value */
/*     ------  ---------  ----------------------------------------------- */
/*     H0      RWORK(5)   Step size to be attempted on the first step. */
/*                        The default value is determined by the solver. */
/*     HMAX    RWORK(6)   Maximum absolute step size allowed.  The */
/*                        default value is infinite. */
/*     HMIN    RWORK(7)   Minimum absolute step size allowed.  The */
/*                        default value is 0.  (This lower bound is not */
/*                        enforced on the final step before reaching */
/*                        TCRIT when ITASK = 4 or 5.) */
/*     MAXORD  IWORK(5)   Maximum order to be allowed.  The default value */
/*                        is 12 if METH = 1, and 5 if METH = 2. (See the */
/*                        MF description above for METH.)  If MAXORD */
/*                        exceeds the default value, it will be reduced */
/*                        to the default value.  If MAXORD is changed */
/*                        during the problem, it may cause the current */
/*                        order to be reduced. */
/*     MXSTEP  IWORK(6)   Maximum number of (internally defined) steps */
/*                        allowed during one call to the solver.  The */
/*                        default value is 500. */
/*     MXHNIL  IWORK(7)   Maximum number of messages printed (per */
/*                        problem) warning that T + H = T on a step */
/*                        (H = step size).  This must be positive to */
/*                        result in a nondefault value.  The default */
/*                        value is 10. */

/*     Optional Outputs */
/*     ---------------- */
/*     As optional additional output from SLSODE, the variables listed */
/*     below are quantities related to the performance of SLSODE which */
/*     are available to the user.  These are communicated by way of the */
/*     work arrays, but also have internal mnemonic names as shown. */
/*     Except where stated otherwise, all of these outputs are defined on */
/*     any successful return from SLSODE, and on any return with ISTATE = */
/*     -1, -2, -4, -5, or -6.  On an illegal input return (ISTATE = -3), */
/*     they will be unchanged from their existing values (if any), except */
/*     possibly for TOLSF, LENRW, and LENIW.  On any error return, */
/*     outputs relevant to the error will be defined, as noted below. */

/*     Name   Location   Meaning */
/*     -----  ---------  ------------------------------------------------ */
/*     HU     RWORK(11)  Step size in t last used (successfully). */
/*     HCUR   RWORK(12)  Step size to be attempted on the next step. */
/*     TCUR   RWORK(13)  Current value of the independent variable which */
/*                       the solver has actually reached, i.e., the */
/*                       current internal mesh point in t. On output, */
/*                       TCUR will always be at least as far as the */
/*                       argument T, but may be farther (if interpolation */
/*                       was done). */
/*     TOLSF  RWORK(14)  Tolerance scale factor, greater than 1.0, */
/*                       computed when a request for too much accuracy */
/*                       was detected (ISTATE = -3 if detected at the */
/*                       start of the problem, ISTATE = -2 otherwise). */
/*                       If ITOL is left unaltered but RTOL and ATOL are */
/*                       uniformly scaled up by a factor of TOLSF for the */
/*                       next call, then the solver is deemed likely to */
/*                       succeed.  (The user may also ignore TOLSF and */
/*                       alter the tolerance parameters in any other way */
/*                       appropriate.) */
/*     NST    IWORK(11)  Number of steps taken for the problem so far. */
/*     NFE    IWORK(12)  Number of F evaluations for the problem so far. */
/*     NJE    IWORK(13)  Number of Jacobian evaluations (and of matrix LU */
/*                       decompositions) for the problem so far. */
/*     NQU    IWORK(14)  Method order last used (successfully). */
/*     NQCUR  IWORK(15)  Order to be attempted on the next step. */
/*     IMXER  IWORK(16)  Index of the component of largest magnitude in */
/*                       the weighted local error vector ( e(i)/EWT(i) ), */
/*                       on an error return with ISTATE = -4 or -5. */
/*     LENRW  IWORK(17)  Length of RWORK actually required.  This is */
/*                       defined on normal returns and on an illegal */
/*                       input return for insufficient storage. */
/*     LENIW  IWORK(18)  Length of IWORK actually required.  This is */
/*                       defined on normal returns and on an illegal */
/*                       input return for insufficient storage. */

/*     The following two arrays are segments of the RWORK array which may */
/*     also be of interest to the user as optional outputs.  For each */
/*     array, the table below gives its internal name, its base address */
/*     in RWORK, and its description. */

/*     Name  Base address  Description */
/*     ----  ------------  ---------------------------------------------- */
/*     YH    21            The Nordsieck history array, of size NYH by */
/*                         (NQCUR + 1), where NYH is the initial value of */
/*                         NEQ.  For j = 0,1,...,NQCUR, column j + 1 of */
/*                         YH contains HCUR**j/factorial(j) times the jth */
/*                         derivative of the interpolating polynomial */
/*                         currently representing the solution, evaluated */
/*                         at t = TCUR. */
/*     ACOR  LENRW-NEQ+1   Array of size NEQ used for the accumulated */
/*                         corrections on each step, scaled on output to */
/*                         represent the estimated local error in Y on */
/*                         the last step.  This is the vector e in the */
/*                         description of the error control.  It is */
/*                         defined only on successful return from SLSODE. */


/*                    Part 2.  Other Callable Routines */
/*                    -------------------------------- */

/*     The following are optional calls which the user may make to gain */
/*     additional capabilities in conjunction with SLSODE. */

/*     Form of call              Function */
/*     ------------------------  ---------------------------------------- */
/*     CALL XSETUN(LUN)          Set the logical unit number, LUN, for */
/*                               output of messages from SLSODE, if the */
/*                               default is not desired.  The default */
/*                               value of LUN is 6. This call may be made */
/*                               at any time and will take effect */
/*                               immediately. */
/*     CALL XSETF(MFLAG)         Set a flag to control the printing of */
/*                               messages by SLSODE.  MFLAG = 0 means do */
/*                               not print.  (Danger:  this risks losing */
/*                               valuable information.)  MFLAG = 1 means */
/*                               print (the default).  This call may be */
/*                               made at any time and will take effect */
/*                               immediately. */
/*     CALL SSRCOM(RSAV,ISAV,JOB)  Saves and restores the contents of the */
/*                               internal COMMON blocks used by SLSODE */
/*                               (see Part 3 below).  RSAV must be a */
/*                               real array of length 218 or more, and */
/*                               ISAV must be an integer array of length */
/*                               37 or more.  JOB = 1 means save COMMON */
/*                               into RSAV/ISAV.  JOB = 2 means restore */
/*                               COMMON from same.  SSRCOM is useful if */
/*                               one is interrupting a run and restarting */
/*                               later, or alternating between two or */
/*                               more problems solved with SLSODE. */
/*     CALL SINTDY(,,,,,)        Provide derivatives of y, of various */
/*     (see below)               orders, at a specified point t, if */
/*                               desired.  It may be called only after a */
/*                               successful return from SLSODE.  Detailed */
/*                               instructions follow. */

/*     Detailed instructions for using SINTDY */
/*     -------------------------------------- */
/*     The form of the CALL is: */

/*           CALL SINTDY (T, K, RWORK(21), NYH, DKY, IFLAG) */

/*     The input parameters are: */

/*     T          Value of independent variable where answers are */
/*                desired (normally the same as the T last returned by */
/*                SLSODE).  For valid results, T must lie between */
/*                TCUR - HU and TCUR.  (See "Optional Outputs" above */
/*                for TCUR and HU.) */
/*     K          Integer order of the derivative desired.  K must */
/*                satisfy 0 <= K <= NQCUR, where NQCUR is the current */
/*                order (see "Optional Outputs").  The capability */
/*                corresponding to K = 0, i.e., computing y(t), is */
/*                already provided by SLSODE directly.  Since */
/*                NQCUR >= 1, the first derivative dy/dt is always */
/*                available with SINTDY. */
/*     RWORK(21)  The base address of the history array YH. */
/*     NYH        Column length of YH, equal to the initial value of NEQ. */

/*     The output parameters are: */

/*     DKY        Real array of length NEQ containing the computed value */
/*                of the Kth derivative of y(t). */
/*     IFLAG      Integer flag, returned as 0 if K and T were legal, */
/*                -1 if K was illegal, and -2 if T was illegal. */
/*                On an error return, a message is also written. */


/*                          Part 3.  Common Blocks */
/*                          ---------------------- */

/*     If SLSODE is to be used in an overlay situation, the user must */
/*     declare, in the primary overlay, the variables in: */
/*     (1) the call sequence to SLSODE, */
/*     (2) the internal COMMON block /SLS001/, of length 255 */
/*         (218 single precision words followed by 37 integer words). */

/*     If SLSODE is used on a system in which the contents of internal */
/*     COMMON blocks are not preserved between calls, the user should */
/*     declare the above COMMON block in his main program to insure that */
/*     its contents are preserved. */

/*     If the solution of a given problem by SLSODE is to be interrupted */
/*     and then later continued, as when restarting an interrupted run or */
/*     alternating between two or more problems, the user should save, */
/*     following the return from the last SLSODE call prior to the */
/*     interruption, the contents of the call sequence variables and the */
/*     internal COMMON block, and later restore these values before the */
/*     next SLSODE call for that problem.   In addition, if XSETUN and/or */
/*     XSETF was called for non-default handling of error messages, then */
/*     these calls must be repeated.  To save and restore the COMMON */
/*     block, use subroutine SSRCOM (see Part 2 above). */


/*              Part 4.  Optionally Replaceable Solver Routines */
/*              ----------------------------------------------- */

/*     Below are descriptions of two routines in the SLSODE package which */
/*     relate to the measurement of errors.  Either routine can be */
/*     replaced by a user-supplied version, if desired.  However, since */
/*     such a replacement may have a major impact on performance, it */
/*     should be done only when absolutely necessary, and only with great */
/*     caution.  (Note:  The means by which the package version of a */
/*     routine is superseded by the user's version may be system- */
/*     dependent.) */

/*     SEWSET */
/*     ------ */
/*     The following subroutine is called just before each internal */
/*     integration step, and sets the array of error weights, EWT, as */
/*     described under ITOL/RTOL/ATOL above: */

/*           SUBROUTINE SEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT) */

/*     where NEQ, ITOL, RTOL, and ATOL are as in the SLSODE call */
/*     sequence, YCUR contains the current dependent variable vector, */
/*     and EWT is the array of weights set by SEWSET. */

/*     If the user supplies this subroutine, it must return in EWT(i) */
/*     (i = 1,...,NEQ) a positive quantity suitable for comparing errors */
/*     in Y(i) to.  The EWT array returned by SEWSET is passed to the */
/*     SVNORM routine (see below), and also used by SLSODE in the */
/*     computation of the optional output IMXER, the diagonal Jacobian */
/*     approximation, and the increments for difference quotient */
/*     Jacobians. */

/*     In the user-supplied version of SEWSET, it may be desirable to use */
/*     the current values of derivatives of y. Derivatives up to order NQ */
/*     are available from the history array YH, described above under */
/*     optional outputs.  In SEWSET, YH is identical to the YCUR array, */
/*     extended to NQ + 1 columns with a column length of NYH and scale */
/*     factors of H**j/factorial(j).  On the first call for the problem, */
/*     given by NST = 0, NQ is 1 and H is temporarily set to 1.0. */
/*     NYH is the initial value of NEQ.  The quantities NQ, H, and NST */
/*     can be obtained by including in SEWSET the statements: */
/*           REAL RLS */
/*           COMMON /SLS001/ RLS(218),ILS(37) */
/*           NQ = ILS(33) */
/*           NST = ILS(34) */
/*           H = RLS(212) */
/*     Thus, for example, the current value of dy/dt can be obtained as */
/*     YCUR(NYH+i)/H (i=1,...,NEQ) (and the division by H is unnecessary */
/*     when NST = 0). */

/*     SVNORM */
/*     ------ */
/*     SVNORM is a real function routine which computes the weighted */
/*     root-mean-square norm of a vector v: */

/*        d = SVNORM (n, v, w) */

/*     where: */
/*     n = the length of the vector, */
/*     v = real array of length n containing the vector, */
/*     w = real array of length n containing weights, */
/*     d = SQRT( (1/n) * sum(v(i)*w(i))**2 ). */

/*     SVNORM is called with n = NEQ and with w(i) = 1.0/EWT(i), where */
/*     EWT is as set by subroutine SEWSET. */

/*     If the user supplies this function, it should return a nonnegative */
/*     value of SVNORM suitable for use in the error control in SLSODE. */
/*     None of the arguments should be altered by SVNORM.  For example, a */
/*     user-supplied SVNORM routine might: */
/*     - Substitute a max-norm of (v(i)*w(i)) for the rms-norm, or */
/*     - Ignore some components of v in the norm, with the effect of */
/*       suppressing the error control on those components of Y. */
/*  --------------------------------------------------------------------- */
/* ***ROUTINES CALLED  SEWSET, SINTDY, RUMACH, SSTODE, SVNORM, XERRWV */
/* ***COMMON BLOCKS    SLS001 */
/* ***REVISION HISTORY  (YYYYMMDD) */
/* 19791129  DATE WRITTEN */
/* 19791213  Minor changes to declarations; DELP init. in STODE. */
/* 19800118  Treat NEQ as array; integer declarations added throughout; */
/*           minor changes to prologue. */
/* 19800306  Corrected TESCO(1,NQP1) setting in CFODE. */
/* 19800519  Corrected access of YH on forced order reduction; */
/*           numerous corrections to prologues and other comments. */
/* 19800617  In main driver, added loading of SQRT(UROUND) in RWORK; */
/*           minor corrections to main prologue. */
/* 19800923  Added zero initialization of HU and NQU. */
/* 19801218  Revised XERRWV routine; minor corrections to main prologue. */
/* 19810401  Minor changes to comments and an error message. */
/* 19810814  Numerous revisions: replaced EWT by 1/EWT; used flags */
/*           JCUR, ICF, IERPJ, IERSL between STODE and subordinates; */
/*           added tuning parameters CCMAX, MAXCOR, MSBP, MXNCF; */
/*           reorganized returns from STODE; reorganized type decls.; */
/*           fixed message length in XERRWV; changed default LUNIT to 6; */
/*           changed Common lengths; changed comments throughout. */
/* 19870330  Major update by ACH: corrected comments throughout; */
/*           removed TRET from Common; rewrote EWSET with 4 loops; */
/*           fixed t test in INTDY; added Cray directives in STODE; */
/*           in STODE, fixed DELP init. and logic around PJAC call; */
/*           combined routines to save/restore Common; */
/*           passed LEVEL = 0 in error message calls (except run abort). */
/* 19890426  Modified prologue to SLATEC/LDOC format.  (FNF) */
/* 19890501  Many improvements to prologue.  (FNF) */
/* 19890503  A few final corrections to prologue.  (FNF) */
/* 19890504  Minor cosmetic changes.  (FNF) */
/* 19890510  Corrected description of Y in Arguments section.  (FNF) */
/* 19890517  Minor corrections to prologue.  (FNF) */
/* 19920514  Updated with prologue edited 891025 by G. Shaw for manual. */
/* 19920515  Converted source lines to upper case.  (FNF) */
/* 19920603  Revised XERRWV calls using mixed upper-lower case.  (ACH) */
/* 19920616  Revised prologue comment regarding CFT.  (ACH) */
/* 19921116  Revised prologue comments regarding Common.  (ACH). */
/* 19930326  Added comment about non-reentrancy.  (FNF) */
/* 19930723  Changed R1MACH to RUMACH. (FNF) */
/* 19930801  Removed ILLIN and NTREP from Common (affects driver logic); */
/*           minor changes to prologue and internal comments; */
/*           changed Hollerith strings to quoted strings; */
/*           changed internal comments to mixed case; */
/*           replaced XERRWV with new version using character type; */
/*           changed dummy dimensions from 1 to *. (ACH) */
/* 19930809  Changed to generic intrinsic names; changed names of */
/*           subprograms and Common blocks to SLSODE etc. (ACH) */
/* 19930929  Eliminated use of REAL intrinsic; other minor changes. (ACH) */
/* 20010412  Removed all 'own' variables from Common block /SLS001/ */
/*           (affects declarations in 6 routines). (ACH) */
/* 20010509  Minor corrections to prologue. (ACH) */
/* 20031105  Restored 'own' variables to Common block /SLS001/, to */
/*           enable interrupt/restart feature. (ACH) */
/* 20031112  Added SAVE statements for data-loaded constants. */

/* ***  END PROLOGUE  SLSODE */

/* *Internal Notes: */

/* Other Routines in the SLSODE Package. */

/* In addition to Subroutine SLSODE, the SLSODE package includes the */
/* following subroutines and function routines: */
/*  SINTDY   computes an interpolated value of the y vector at t = TOUT. */
/*  SSTODE   is the core integrator, which does one step of the */
/*           integration and the associated error control. */
/*  SCFODE   sets all method coefficients and test constants. */
/*  SPREPJ   computes and preprocesses the Jacobian matrix J = df/dy */
/*           and the Newton iteration matrix P = I - h*l0*J. */
/*  SSOLSY   manages solution of linear system in chord iteration. */
/*  SEWSET   sets the error weight vector EWT before each step. */
/*  SVNORM   computes the weighted R.M.S. norm of a vector. */
/*  SSRCOM   is a user-callable routine to save and restore */
/*           the contents of the internal Common block. */
/*  SGEFA and SGESL   are routines from LINPACK for solving full */
/*           systems of linear algebraic equations. */
/*  SGBFA and SGBSL   are routines from LINPACK for solving banded */
/*           linear systems. */
/*  RUMACH   computes the unit roundoff in a machine-independent manner. */
/*  XERRWV, XSETUN, XSETF, IXSAV, IUMACH   handle the printing of all */
/*           error messages and warnings.  XERRWV is machine-dependent. */
/* Note: SVNORM, RUMACH, IXSAV, and IUMACH are function routines. */
/* All the others are subroutines. */

/* **End */

/*  Declare externals. */

/*  Declare all other variables. */
/* ----------------------------------------------------------------------- */
/* The following internal Common block contains */
/* (a) variables which are local to any subroutine but whose values must */
/*     be preserved between calls to the routine ("own" variables), and */
/* (b) variables which are communicated between subroutines. */
/* The block SLS001 is declared in subroutines SLSODE, SINTDY, SSTODE, */
/* SPREPJ, and SSOLSY. */
/* Groups of variables are replaced by dummy arrays in the Common */
/* declarations in routines where those variables are not used. */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --neq;
    --y;
    --rtol;
    --atol;
    --rwork;
    --iwork;

    /* Function Body */
/* ----------------------------------------------------------------------- */
/* Block A. */
/* This code block is executed on every call. */
/* It tests ISTATE and ITASK for legality and branches appropriately. */
/* If ISTATE .GT. 1 but the flag INIT shows that initialization has */
/* not yet been done, an error return occurs. */
/* If ISTATE = 1 and TOUT = T, return immediately. */
/* ----------------------------------------------------------------------- */

/* ***FIRST EXECUTABLE STATEMENT  SLSODE */
    if (*istate < 1 || *istate > 3) {
	goto L601;
    }
    if (*itask < 1 || *itask > 5) {
	goto L602;
    }
    if (*istate == 1) {
	goto L10;
    }
    if (sls001_1.init == 0) {
	goto L603;
    }
    if (*istate == 2) {
	goto L200;
    }
    goto L20;
L10:
    sls001_1.init = 0;
    if (*tout == *t) {
	return 0;
    }
/* ----------------------------------------------------------------------- */
/* Block B. */
/* The next code block is executed for the initial call (ISTATE = 1), */
/* or for a continuation call with parameter changes (ISTATE = 3). */
/* It contains checking of all inputs and various initializations. */

/* First check legality of the non-optional inputs NEQ, ITOL, IOPT, */
/* MF, ML, and MU. */
/* ----------------------------------------------------------------------- */
L20:
    if (neq[1] <= 0) {
	goto L604;
    }
    if (*istate == 1) {
	goto L25;
    }
    if (neq[1] > sls001_1.n) {
	goto L605;
    }
L25:
    sls001_1.n = neq[1];
    if (*itol < 1 || *itol > 4) {
	goto L606;
    }
    if (*iopt < 0 || *iopt > 1) {
	goto L607;
    }
    sls001_1.meth = *mf / 10;
    sls001_1.miter = *mf - sls001_1.meth * 10;
    if (sls001_1.meth < 1 || sls001_1.meth > 2) {
	goto L608;
    }
    if (sls001_1.miter < 0 || sls001_1.miter > 5) {
	goto L608;
    }
    if (sls001_1.miter <= 3) {
	goto L30;
    }
    ml = iwork[1];
    mu = iwork[2];
    if (ml < 0 || ml >= sls001_1.n) {
	goto L609;
    }
    if (mu < 0 || mu >= sls001_1.n) {
	goto L610;
    }
L30:
/* Next process and check the optional inputs. -------------------------- */
    if (*iopt == 1) {
	goto L40;
    }
    sls001_1.maxord = mord[sls001_1.meth - 1];
    sls001_1.mxstep = mxstp0;
    sls001_1.mxhnil = mxhnl0;
    if (*istate == 1) {
	h0 = 0.f;
    }
    sls001_1.hmxi = 0.f;
    sls001_1.hmin = 0.f;
    goto L60;
L40:
    sls001_1.maxord = iwork[5];
    if (sls001_1.maxord < 0) {
	goto L611;
    }
    if (sls001_1.maxord == 0) {
	sls001_1.maxord = 100;
    }
/* Computing MIN */
    i__1 = sls001_1.maxord, i__2 = mord[sls001_1.meth - 1];
    sls001_1.maxord = min(i__1,i__2);
    sls001_1.mxstep = iwork[6];
    if (sls001_1.mxstep < 0) {
	goto L612;
    }
    if (sls001_1.mxstep == 0) {
	sls001_1.mxstep = mxstp0;
    }
    sls001_1.mxhnil = iwork[7];
    if (sls001_1.mxhnil < 0) {
	goto L613;
    }
    if (sls001_1.mxhnil == 0) {
	sls001_1.mxhnil = mxhnl0;
    }
    if (*istate != 1) {
	goto L50;
    }
    h0 = rwork[5];
    if ((*tout - *t) * h0 < 0.f) {
	goto L614;
    }
L50:
    hmax = rwork[6];
    if (hmax < 0.f) {
	goto L615;
    }
    sls001_1.hmxi = 0.f;
    if (hmax > 0.f) {
	sls001_1.hmxi = 1.f / hmax;
    }
    sls001_1.hmin = rwork[7];
    if (sls001_1.hmin < 0.f) {
	goto L616;
    }
/* ----------------------------------------------------------------------- */
/* Set work array pointers and check lengths LRW and LIW. */
/* Pointers to segments of RWORK and IWORK are named by prefixing L to */
/* the name of the segment.  E.g., the segment YH starts at RWORK(LYH). */
/* Segments of RWORK (in order) are denoted  YH, WM, EWT, SAVF, ACOR. */
/* ----------------------------------------------------------------------- */
L60:
    sls001_1.lyh = 21;
    if (*istate == 1) {
	sls001_1.nyh = sls001_1.n;
    }
    sls001_1.lwm = sls001_1.lyh + (sls001_1.maxord + 1) * sls001_1.nyh;
    if (sls001_1.miter == 0) {
	lenwm = 0;
    }
    if (sls001_1.miter == 1 || sls001_1.miter == 2) {
	lenwm = sls001_1.n * sls001_1.n + 2;
    }
    if (sls001_1.miter == 3) {
	lenwm = sls001_1.n + 2;
    }
    if (sls001_1.miter >= 4) {
	lenwm = ((ml << 1) + mu + 1) * sls001_1.n + 2;
    }
    sls001_1.lewt = sls001_1.lwm + lenwm;
    sls001_1.lsavf = sls001_1.lewt + sls001_1.n;
    sls001_1.lacor = sls001_1.lsavf + sls001_1.n;
    lenrw = sls001_1.lacor + sls001_1.n - 1;
    iwork[17] = lenrw;
    sls001_1.liwm = 1;
    leniw = sls001_1.n + 20;
    if (sls001_1.miter == 0 || sls001_1.miter == 3) {
	leniw = 20;
    }
    iwork[18] = leniw;
    if (lenrw > *lrw) {
	goto L617;
    }
    if (leniw > *liw) {
	goto L618;
    }
/* Check RTOL and ATOL for legality. ------------------------------------ */
    rtoli = rtol[1];
    atoli = atol[1];
    i__1 = sls001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (*itol >= 3) {
	    rtoli = rtol[i__];
	}
	if (*itol == 2 || *itol == 4) {
	    atoli = atol[i__];
	}
	if (rtoli < 0.f) {
	    goto L619;
	}
	if (atoli < 0.f) {
	    goto L620;
	}
/* L70: */
    }
    if (*istate == 1) {
	goto L100;
    }
/* If ISTATE = 3, set flag to signal parameter changes to SSTODE. ------- */
    sls001_1.jstart = -1;
    if (sls001_1.nq <= sls001_1.maxord) {
	goto L90;
    }
/* MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into SAVF. --------- */
    i__1 = sls001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L80: */
	rwork[i__ + sls001_1.lsavf - 1] = rwork[i__ + sls001_1.lwm - 1];
    }
/* Reload WM(1) = RWORK(LWM), since LWM may have changed. --------------- */
L90:
    if (sls001_1.miter > 0) {
	rwork[sls001_1.lwm] = sqrt(sls001_1.uround);
    }
    if (sls001_1.n == sls001_1.nyh) {
	goto L200;
    }
/* NEQ was reduced.  Zero part of YH to avoid undefined references. ----- */
    i1 = sls001_1.lyh + sls001_1.l * sls001_1.nyh;
    i2 = sls001_1.lyh + (sls001_1.maxord + 1) * sls001_1.nyh - 1;
    if (i1 > i2) {
	goto L200;
    }
    i__1 = i2;
    for (i__ = i1; i__ <= i__1; ++i__) {
/* L95: */
	rwork[i__] = 0.f;
    }
    goto L200;
/* ----------------------------------------------------------------------- */
/* Block C. */
/* The next block is for the initial call only (ISTATE = 1). */
/* It contains all remaining initializations, the initial call to F, */
/* and the calculation of the initial step size. */
/* The error weights in EWT are inverted after being loaded. */
/* ----------------------------------------------------------------------- */
L100:
    sls001_1.uround = rumach_();
    sls001_1.tn = *t;
    if (*itask != 4 && *itask != 5) {
	goto L110;
    }
    tcrit = rwork[1];
    if ((tcrit - *tout) * (*tout - *t) < 0.f) {
	goto L625;
    }
    if (h0 != 0.f && (*t + h0 - tcrit) * h0 > 0.f) {
	h0 = tcrit - *t;
    }
L110:
    sls001_1.jstart = 0;
    if (sls001_1.miter > 0) {
	rwork[sls001_1.lwm] = sqrt(sls001_1.uround);
    }
    sls001_1.nhnil = 0;
    sls001_1.nst = 0;
    sls001_1.nje = 0;
    sls001_1.nslast = 0;
    sls001_1.hu = 0.f;
    sls001_1.nqu = 0;
    sls001_1.ccmax = .3f;
    sls001_1.maxcor = 3;
    sls001_1.msbp = 20;
    sls001_1.mxncf = 10;
/* Initial call to F.  (LF0 points to YH(*,2).) ------------------------- */
    lf0 = sls001_1.lyh + sls001_1.nyh;
    (*f)(&neq[1], t, &y[1], &rwork[lf0]);
    sls001_1.nfe = 1;
/* Load the initial value vector in YH. --------------------------------- */
    i__1 = sls001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L115: */
	rwork[i__ + sls001_1.lyh - 1] = y[i__];
    }
/* Load and invert the EWT array.  (H is temporarily set to 1.0.) ------- */
    sls001_1.nq = 1;
    sls001_1.h__ = 1.f;
    sewset_(&sls001_1.n, itol, &rtol[1], &atol[1], &rwork[sls001_1.lyh], &
	    rwork[sls001_1.lewt]);
    i__1 = sls001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (rwork[i__ + sls001_1.lewt - 1] <= 0.f) {
	    goto L621;
	}
/* L120: */
	rwork[i__ + sls001_1.lewt - 1] = 1.f / rwork[i__ + sls001_1.lewt - 1];
    }
/* ----------------------------------------------------------------------- */
/* The coding below computes the step size, H0, to be attempted on the */
/* first step, unless the user has supplied a value for this. */
/* First check that TOUT - T differs significantly from zero. */
/* A scalar tolerance quantity TOL is computed, as MAX(RTOL(I)) */
/* if this is positive, or MAX(ATOL(I)/ABS(Y(I))) otherwise, adjusted */
/* so as to be between 100*UROUND and 1.0E-3. */
/* Then the computed value H0 is given by.. */
/*                                      NEQ */
/*   H0**2 = TOL / ( w0**-2 + (1/NEQ) * SUM ( f(i)/ywt(i) )**2  ) */
/*                                       1 */
/* where   w0     = MAX ( ABS(T), ABS(TOUT) ), */
/*         f(i)   = i-th component of initial value of f, */
/*         ywt(i) = EWT(i)/TOL  (a weight for y(i)). */
/* The sign of H0 is inferred from the initial values of TOUT and T. */
/* ----------------------------------------------------------------------- */
    if (h0 != 0.f) {
	goto L180;
    }
    tdist = (r__1 = *tout - *t, dabs(r__1));
/* Computing MAX */
    r__1 = dabs(*t), r__2 = dabs(*tout);
    w0 = dmax(r__1,r__2);
    if (tdist < sls001_1.uround * 2.f * w0) {
	goto L622;
    }
    tol = rtol[1];
    if (*itol <= 2) {
	goto L140;
    }
    i__1 = sls001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L130: */
/* Computing MAX */
	r__1 = tol, r__2 = rtol[i__];
	tol = dmax(r__1,r__2);
    }
L140:
    if (tol > 0.f) {
	goto L160;
    }
    atoli = atol[1];
    i__1 = sls001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (*itol == 2 || *itol == 4) {
	    atoli = atol[i__];
	}
	ayi = (r__1 = y[i__], dabs(r__1));
	if (ayi != 0.f) {
/* Computing MAX */
	    r__1 = tol, r__2 = atoli / ayi;
	    tol = dmax(r__1,r__2);
	}
/* L150: */
    }
L160:
/* Computing MAX */
    r__1 = tol, r__2 = sls001_1.uround * 100.f;
    tol = dmax(r__1,r__2);
    tol = dmin(tol,.001f);
    sum = svnorm_(&sls001_1.n, &rwork[lf0], &rwork[sls001_1.lewt]);
/* Computing 2nd power */
    r__1 = sum;
    sum = 1.f / (tol * w0 * w0) + tol * (r__1 * r__1);
    h0 = 1.f / sqrt(sum);
    h0 = dmin(h0,tdist);
    r__1 = *tout - *t;
    h0 = r_sign(&h0, &r__1);
/* Adjust H0 if necessary to meet HMAX bound. --------------------------- */
L180:
    rh = dabs(h0) * sls001_1.hmxi;
    if (rh > 1.f) {
	h0 /= rh;
    }
/* Load H with H0 and scale YH(*,2) by H0. ------------------------------ */
    sls001_1.h__ = h0;
    i__1 = sls001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L190: */
	rwork[i__ + lf0 - 1] = h0 * rwork[i__ + lf0 - 1];
    }
    goto L270;
/* ----------------------------------------------------------------------- */
/* Block D. */
/* The next code block is for continuation calls only (ISTATE = 2 or 3) */
/* and is to check stop conditions before taking a step. */
/* ----------------------------------------------------------------------- */
L200:
    sls001_1.nslast = sls001_1.nst;
    switch (*itask) {
	case 1:  goto L210;
	case 2:  goto L250;
	case 3:  goto L220;
	case 4:  goto L230;
	case 5:  goto L240;
    }
L210:
    if ((sls001_1.tn - *tout) * sls001_1.h__ < 0.f) {
	goto L250;
    }
    sintdy_(tout, &c__0, &rwork[sls001_1.lyh], &sls001_1.nyh, &y[1], &iflag);
    if (iflag != 0) {
	goto L627;
    }
    *t = *tout;
    goto L420;
L220:
    tp = sls001_1.tn - sls001_1.hu * (sls001_1.uround * 100.f + 1.f);
    if ((tp - *tout) * sls001_1.h__ > 0.f) {
	goto L623;
    }
    if ((sls001_1.tn - *tout) * sls001_1.h__ < 0.f) {
	goto L250;
    }
    goto L400;
L230:
    tcrit = rwork[1];
    if ((sls001_1.tn - tcrit) * sls001_1.h__ > 0.f) {
	goto L624;
    }
    if ((tcrit - *tout) * sls001_1.h__ < 0.f) {
	goto L625;
    }
    if ((sls001_1.tn - *tout) * sls001_1.h__ < 0.f) {
	goto L245;
    }
    sintdy_(tout, &c__0, &rwork[sls001_1.lyh], &sls001_1.nyh, &y[1], &iflag);
    if (iflag != 0) {
	goto L627;
    }
    *t = *tout;
    goto L420;
L240:
    tcrit = rwork[1];
    if ((sls001_1.tn - tcrit) * sls001_1.h__ > 0.f) {
	goto L624;
    }
L245:
    hmx = dabs(sls001_1.tn) + dabs(sls001_1.h__);
    ihit = (r__1 = sls001_1.tn - tcrit, dabs(r__1)) <= sls001_1.uround * 
	    100.f * hmx;
    if (ihit) {
	goto L400;
    }
    tnext = sls001_1.tn + sls001_1.h__ * (sls001_1.uround * 4.f + 1.f);
    if ((tnext - tcrit) * sls001_1.h__ <= 0.f) {
	goto L250;
    }
    sls001_1.h__ = (tcrit - sls001_1.tn) * (1.f - sls001_1.uround * 4.f);
    if (*istate == 2) {
	sls001_1.jstart = -2;
    }
/* ----------------------------------------------------------------------- */
/* Block E. */
/* The next block is normally executed for all calls and contains */
/* the call to the one-step core integrator SSTODE. */

/* This is a looping point for the integration steps. */

/* First check for too many steps being taken, update EWT (if not at */
/* start of problem), check for too much accuracy being requested, and */
/* check for H below the roundoff level in T. */
/* ----------------------------------------------------------------------- */
L250:
    if (sls001_1.nst - sls001_1.nslast >= sls001_1.mxstep) {
	goto L500;
    }
    sewset_(&sls001_1.n, itol, &rtol[1], &atol[1], &rwork[sls001_1.lyh], &
	    rwork[sls001_1.lewt]);
    i__1 = sls001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (rwork[i__ + sls001_1.lewt - 1] <= 0.f) {
	    goto L510;
	}
/* L260: */
	rwork[i__ + sls001_1.lewt - 1] = 1.f / rwork[i__ + sls001_1.lewt - 1];
    }
L270:
    tolsf = sls001_1.uround * svnorm_(&sls001_1.n, &rwork[sls001_1.lyh], &
	    rwork[sls001_1.lewt]);
    if (tolsf <= 1.f) {
	goto L280;
    }
    tolsf *= 2.f;
    if (sls001_1.nst == 0) {
	goto L626;
    }
    goto L520;
L280:
    if (sls001_1.tn + sls001_1.h__ != sls001_1.tn) {
	goto L290;
    }
    ++sls001_1.nhnil;
    if (sls001_1.nhnil > sls001_1.mxhnil) {
	goto L290;
    }
    if (sls001_1.nhnil < sls001_1.mxhnil) {
	goto L290;
    }
L290:
/* ----------------------------------------------------------------------- */
/*  CALL SSTODE(NEQ,Y,YH,NYH,YH,EWT,SAVF,ACOR,WM,IWM,F,JAC,SPREPJ,SSOLSY) */
/* ----------------------------------------------------------------------- */
    sstode_(&neq[1], &y[1], &rwork[sls001_1.lyh], &sls001_1.nyh, &rwork[
	    sls001_1.lyh], &rwork[sls001_1.lewt], &rwork[sls001_1.lsavf], &
	    rwork[sls001_1.lacor], &rwork[sls001_1.lwm], &iwork[sls001_1.liwm]
	    , (S_fp)f, (U_fp)jac, (U_fp)sprepj_, (U_fp)ssolsy_);
    kgo = 1 - sls001_1.kflag;
    switch (kgo) {
	case 1:  goto L300;
	case 2:  goto L530;
	case 3:  goto L540;
    }
/* ----------------------------------------------------------------------- */
/* Block F. */
/* The following block handles the case of a successful return from the */
/* core integrator (KFLAG = 0).  Test for stop conditions. */
/* ----------------------------------------------------------------------- */
L300:
    sls001_1.init = 1;
    switch (*itask) {
	case 1:  goto L310;
	case 2:  goto L400;
	case 3:  goto L330;
	case 4:  goto L340;
	case 5:  goto L350;
    }
/* ITASK = 1.  If TOUT has been reached, interpolate. ------------------- */
L310:
    if ((sls001_1.tn - *tout) * sls001_1.h__ < 0.f) {
	goto L250;
    }
    sintdy_(tout, &c__0, &rwork[sls001_1.lyh], &sls001_1.nyh, &y[1], &iflag);
    *t = *tout;
    goto L420;
/* ITASK = 3.  Jump to exit if TOUT was reached. ------------------------ */
L330:
    if ((sls001_1.tn - *tout) * sls001_1.h__ >= 0.f) {
	goto L400;
    }
    goto L250;
/* ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary. */
L340:
    if ((sls001_1.tn - *tout) * sls001_1.h__ < 0.f) {
	goto L345;
    }
    sintdy_(tout, &c__0, &rwork[sls001_1.lyh], &sls001_1.nyh, &y[1], &iflag);
    *t = *tout;
    goto L420;
L345:
    hmx = dabs(sls001_1.tn) + dabs(sls001_1.h__);
    ihit = (r__1 = sls001_1.tn - tcrit, dabs(r__1)) <= sls001_1.uround * 
	    100.f * hmx;
    if (ihit) {
	goto L400;
    }
    tnext = sls001_1.tn + sls001_1.h__ * (sls001_1.uround * 4.f + 1.f);
    if ((tnext - tcrit) * sls001_1.h__ <= 0.f) {
	goto L250;
    }
    sls001_1.h__ = (tcrit - sls001_1.tn) * (1.f - sls001_1.uround * 4.f);
    sls001_1.jstart = -2;
    goto L250;
/* ITASK = 5.  See if TCRIT was reached and jump to exit. --------------- */
L350:
    hmx = dabs(sls001_1.tn) + dabs(sls001_1.h__);
    ihit = (r__1 = sls001_1.tn - tcrit, dabs(r__1)) <= sls001_1.uround * 
	    100.f * hmx;
/* ----------------------------------------------------------------------- */
/* Block G. */
/* The following block handles all successful returns from SLSODE. */
/* If ITASK .NE. 1, Y is loaded from YH and T is set accordingly. */
/* ISTATE is set to 2, and the optional outputs are loaded into the */
/* work arrays before returning. */
/* ----------------------------------------------------------------------- */
L400:
    i__1 = sls001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L410: */
	y[i__] = rwork[i__ + sls001_1.lyh - 1];
    }
    *t = sls001_1.tn;
    if (*itask != 4 && *itask != 5) {
	goto L420;
    }
    if (ihit) {
	*t = tcrit;
    }
L420:
    *istate = 2;
    rwork[11] = sls001_1.hu;
    rwork[12] = sls001_1.h__;
    rwork[13] = sls001_1.tn;
    iwork[11] = sls001_1.nst;
    iwork[12] = sls001_1.nfe;
    iwork[13] = sls001_1.nje;
    iwork[14] = sls001_1.nqu;
    iwork[15] = sls001_1.nq;
    return 0;
/* ----------------------------------------------------------------------- */
/* Block H. */
/* The following block handles all unsuccessful returns other than */
/* those for illegal input.  First the error message routine is called. */
/* If there was an error test or convergence test failure, IMXER is set. */
/* Then Y is loaded from YH and T is set to TN.  The optional outputs */
/* are loaded into the work arrays before returning. */
/* ----------------------------------------------------------------------- */
/* The maximum number of steps was taken before reaching TOUT. ---------- */
L500:
    *istate = -1;
    goto L580;
/* EWT(I) .LE. 0.0 for some I (not at start of problem). ---------------- */
L510:
    ewti = rwork[sls001_1.lewt + i__ - 1];
    *istate = -6;
    goto L580;
/* Too much accuracy requested for machine precision. ------------------- */
L520:
    rwork[14] = tolsf;
    *istate = -2;
    goto L580;
/* KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. ----- */
L530:
    *istate = -4;
    goto L560;
/* KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ---- */
L540:
    *istate = -5;
/* Compute IMXER if relevant. ------------------------------------------- */
L560:
    big = 0.f;
    imxer = 1;
    i__1 = sls001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	size = (r__1 = rwork[i__ + sls001_1.lacor - 1] * rwork[i__ + 
		sls001_1.lewt - 1], dabs(r__1));
	if (big >= size) {
	    goto L570;
	}
	big = size;
	imxer = i__;
L570:
	;
    }
    iwork[16] = imxer;
/* Set Y vector, T, and optional outputs. ------------------------------- */
L580:
    i__1 = sls001_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L590: */
	y[i__] = rwork[i__ + sls001_1.lyh - 1];
    }
    *t = sls001_1.tn;
    rwork[11] = sls001_1.hu;
    rwork[12] = sls001_1.h__;
    rwork[13] = sls001_1.tn;
    iwork[11] = sls001_1.nst;
    iwork[12] = sls001_1.nfe;
    iwork[13] = sls001_1.nje;
    iwork[14] = sls001_1.nqu;
    iwork[15] = sls001_1.nq;
    return 0;
/* ----------------------------------------------------------------------- */
/* Block I. */
/* The following block handles all error returns due to illegal input */
/* (ISTATE = -3), as detected before calling the core integrator. */
/* First the error message routine is called.  If the illegal input */
/* is a negative ISTATE, the run is aborted (apparent infinite loop). */
/* ----------------------------------------------------------------------- */
L601:
    if (*istate < 0) {
	goto L800;
    }
    goto L700;
L602:
    goto L700;
L603:
    goto L700;
L604:
    goto L700;
L605:
    goto L700;
L606:
    goto L700;
L607:
    goto L700;
L608:
    goto L700;
L609:
    goto L700;
L610:
    goto L700;
L611:
    goto L700;
L612:
    goto L700;
L613:
    goto L700;
L614:
    goto L700;
L615:
    goto L700;
L616:
    goto L700;
L617:
    goto L700;
L618:
    goto L700;
L619:
    goto L700;
L620:
    goto L700;
L621:
    ewti = rwork[sls001_1.lewt + i__ - 1];
    goto L700;
L622:
    goto L700;
L623:
    goto L700;
L624:
    goto L700;
L625:
    goto L700;
L626:
    rwork[14] = tolsf;
    goto L700;
L627:

L700:
    *istate = -3;
    return 0;

L800:
    return 0;
/* ----------------------- END OF SUBROUTINE SLSODE ---------------------- */
} /* slsode_ */

