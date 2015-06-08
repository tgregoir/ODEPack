#ifdef __cplusplus
extern "C" {
#endif

#include "f2c.h"

/* Subroutine */
__device__ int rumsum_(real *a, real *b, real *c__)
{
/*     Routine to force normal storing of A + B, for RUMACH. */
    *c__ = *a + *b;
    return 0;
} /* rumsum_ */

#ifdef __cplusplus
}
#endif

