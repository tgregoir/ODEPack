#include "odesolve.h"

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
	float y[3] = { 1.0f, 0.0f, 0.0f };
	float tout = 0.4f;
	float rtol = 1e-4f;
	float atol[3] = { 1e-6f, 1e-10f, 1e-6f };
	int iout;

	for (iout = 1; iout <= 12; iout++) {
		odesolve(f, 3, y, 0.0f, tout, 2, &rtol, atol, jac, 1);
		printf(" At t=%.4e   y=%.6e, %.6e, %.6e\n",
		       tout, y[0], y[1], y[2]);
		tout = tout * 10.0f;
	}

	return 0;
}
