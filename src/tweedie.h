#ifndef CPLM_TWEEDIE_H 
#define CPLM_TWEEDIE_H 

void dtweedie(int n,  double *y, double *mu, double phi, double p,
              double *wts, double *ans) ;

double dl2tweedie(int n, double *y, double *mu, double phi, double p,
                  double *wts) ;
SEXP cplm_dltweedie(SEXP y, SEXP mu, SEXP phi, SEXP p, SEXP wts) ;

#endif
