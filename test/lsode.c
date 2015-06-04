#include "lsode.h"

#include <stdlib.h>

int odesolve(fct f, int neq, float *y, float ti, float tf,
             int itol, float *rtol, float *atol,
             jacfct jac, int stiff)
{
	/* setting parameters up */
	int itask = 1;
	int istate = 1;
	int iopt = 0;
	int mf = (stiff ? 21 : 10);

	/* compute how much memory we need */
	int lrw;
	int liw;
	if (stiff) {
		lrw = 22 + 9 * neq + neq * neq;
		liw = 20 + neq;
	}
	else {
		lrw = 20 + 16 * neq;
		liw = 20;
	}

	/* memory allocation */
	float *rwork = (float *)malloc(sizeof(float) * lrw);
	if (!rwork) {
		return 1;
	}
	int *iwork = (int *)malloc(sizeof(int) * liw);
	if (!iwork) {
		free(rwork);
		return 1;
	}

	/* solve the ODE */
	slsode_(f, &neq, y, &ti, &tf, &itol, rtol, atol,
	        &itask, &istate, &iopt,
	        rwork, &lrw, iwork, &liw, jac, &mf);

	/* clean up */
	free(iwork);
	free(rwork);

	return (istate < 0);
}
