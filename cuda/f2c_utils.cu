#ifdef __cplusplus
extern "C" {
#endif

#include "f2c.h"
#include <math.h>

__device__ double pow_dd(doublereal *ap, doublereal *bp)
{
	return pow(*ap, *bp);
}

__device__ double pow_ri(real *ap, integer *bp)
{
	double pow, x;
	integer n;
	unsigned long u;

	pow = 1;
	x = *ap;
	n = *bp;

	if (n != 0) {
		if (n < 0) {
			n = -n;
			x = 1 / x;
		}

		for (u = n; ; ) {
			if (u & 01)
				pow *= x;
			if (u >>= 1)
				x *= x;
			else
				break;
		}
	}

	return(pow);
}

__device__ double r_sign(real *a, real *b)
{
	double x;
	x = (*a >= 0 ? *a : -*a);
	return (*b >= 0 ? x : -x);
}

#ifdef __cplusplus
}
#endif

