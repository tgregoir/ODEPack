#ifndef LSODE_H
#define LSODE_H

typedef void (*fct)(int *, float *, float *, float *);
typedef void (*jacfct)(int *, float *, float *, int *, int *, float *,
                       int *);

int slsode_(fct f, int *neq, float *y, float *t, float *tout,
            int *itol, float *rtol, float *atol,
            int *itask, int * istate, int *iopt,
            float *rwork, int *lrw, int *iwork, int *liw,
            jacfct jac, int *mf);

int odesolve(fct, int, float *, float, float, int, float *, float *,
             jacfct, int);

#endif /* LSODE_H */
