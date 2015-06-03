#include "lsode.h"

#include <stdio.h>

void f(int *neq, float *t, float *y, float *ydot)
{
	ydot[0] = -0.04f * y[0] + 1e4f * y[1] * y[2];
	ydot[2] = 3e7f * y[1] * y[1];
	ydot[1] = - ydot[0] - ydot[2];
}

void jac(int *neq, float *t, float *y, int *ml, int *mu,
         float *pd, int *nrpd)
{
        pd[0] = -.04f;
        pd[1] = .04f;
	pd[2] = 0.0f;
        pd[3] = 1e4f * y[2];
	pd[5] = 6e7f * y[1];
	pd[4] = - pd[3] - pd[5];
	pd[6] = 1e4f * y[1];
	pd[7] = -pd[6];
	pd[8] = 0.0f;
}

int main(void)
{
	int neq = 3;
	float y[3] = { 1.0f, 0.0f, 0.0f };
	float t = 0.0f;
	float tout = 0.4f;
	int itol = 2;
	float rtol = 1e-4f;
	float atol[3] = { 1e-6f, 1e-10f, 1e-6f };
	int itask = 1;
	int istate = 1;
	int iopt = 0;
	int lrw = 58;
	int liw = 23;
	int mf = 21;
	float rwork[58];
	int iwork[23];
	int iout;

	for (iout = 1; iout <= 12; iout++) {
		slsode_(f, &neq, y, &t, &tout, &itol, &rtol,
		        atol, &itask, &istate, &iopt, rwork, &lrw,
		        iwork, &liw, jac, &mf);
		printf(" At t=%.4e   y=%.6e, %.6e, %.6e\n",
		       t, y[0], y[1], y[2]);
		if (istate < 0) {
			printf("Error halt.. ISTATE =%d\n", istate);
			return 1;
		}
		tout = tout * 10.0f;
	}

	printf(" No. steps=%d,  No. f-s=%d,  No. J-s=%d\n",
	       iwork[10], iwork[11], iwork[12]);

	return 0;
}
