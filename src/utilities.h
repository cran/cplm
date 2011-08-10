//  
//  Header for utility functions used in cplm
//       Wayne Zhang
//

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <R_ext/Applic.h>

#define EPS 0.001 
#define SCALE (-1) 

double * vect(int n) ;
double ** matrix(int nr, int nc) ;
double cumsum(double *x, int n) ;
double cumwsum(double *x, double *w, int n) ;
void lbfgsbU(double *x, double lower, double upper, double *val,
             optimfn fn, optimgr gr, void *ex,  int *conv) ;
SEXP getListElement (SEXP list, char *str);



