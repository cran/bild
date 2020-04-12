#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(mblik1)(double *logL, double *pij, double *beta, double *lpsi, int *npar, double *x, int *y, double *theta, double *work, int *n);
extern void F77_NAME(mbgd1)(double *gbeta, double *glpsi, double *beta, double *lpsi, int *npar, double *x, int *y, double *theta, double *work, double *der, double *db, double *dbeta, double *dbeta1, int *n);
extern void F77_NAME(integ)(double *logL, double *bt2, double *beta2, double *lpsi, double *omega, int *npar, double *x2, int *y2, double *theta2, double *work2, int *n, double *li, double *ls, double *epsabs, double *epsrel, int *key, int *limit);
extern void F77_NAME(integ1)(double *logL, double *bt2, double *beta2, double *lpsi, double *omega, int *npar, double *x2, int *y2, double *theta2, double *work2, int *n, double *li, double *ls, double *epsabs, double *epsrel, int *key, int *limit);
extern void F77_NAME(gint)(double *gbeta, double *glpsi1, double *glpsi2, double *gvar, double *x2, double *theta2, double *work2, int *y2, double *lpsi,  double *beta2, double *bt2, double *dbt, double *dbt1, double *dbt2, double *dder, double *ddb, double *ddb1, double *ddb2, double *omega, int *npar, int *n, double *li, double *ls, double *epsabs, double *epsrel, int *key, int *limit);
extern void F77_NAME(gint1)(double *gbeta, double *glpsi1, double *gvar, double *bt2, double *beta2, double *lpsi, double *omega, int *npar, double *x2, int *y2, double *theta2, double *work2, int *n, double *dbt, double *dbt1, double *dder, double *ddb, double *li, double *ls, double *epsabs, double *epsrel, int *key, int *limit);
extern void F77_NAME(blik2m)(double *logL, double *pij, double *beta, double *lpsi, int *npar, double *x, int *y, double *theta, double *work, int *n);
extern void F77_NAME(bgd2m)(double *gbeta, double *glpsi1, double *glpsi2, double *beta, double *lpsi, int *npar, double *x, int *y, double *theta, double *work, double *dbeta, double *dbeta1, double *dbeta2, int *n, double *der, double *db, double *db1, double *db2);


static const R_FortranMethodDef FortranEntries[] = {
    {"mblik1",   (DL_FUNC) &F77_NAME(mblik1),   10},
    {"mbgd1",    (DL_FUNC) &F77_NAME(mbgd1),    14},
    {"integ",    (DL_FUNC) &F77_NAME(integ),    17},
    {"integ1",   (DL_FUNC) &F77_NAME(integ1),   17},
    {"gint",     (DL_FUNC) &F77_NAME(gint),     27},
    {"gint1",    (DL_FUNC) &F77_NAME(gint1),    23},
    {"blik2m",   (DL_FUNC) &F77_NAME(blik2m),   10},
    {"bgd2m",    (DL_FUNC) &F77_NAME(bgd2m),    18},

    {NULL, NULL, 0}
};

void R_init_bild(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
