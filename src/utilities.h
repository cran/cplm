//  
//  Header for utility functions used in cplm
//       Wayne Zhang
//

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <R_ext/Applic.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>



#define EPS (0.001) 
#define SCALE (-1) 

double * dvect(int n) ;
double ** dmatrix(int nr, int nc) ;
int * ivect(int n) ;
int ** imatrix(int nr, int nc) ;


double dcumsum(double *x, int n) ;
int icumsum(int *x, int n) ;
double dcumwsum(double *x, double *w, int n) ;
double icumwsum(int *x, double *w, int n) ;

void lbfgsbU(double *x, double lower, double upper, double *val,
             optimfn fn, optimgr gr, void *ex,  int *conv) ;
SEXP getListElement (SEXP list, char *str);

double norm (double *x, int n);
double dist (double *x, double *y, int n);
double varFun (double mu, double p);
double linkFun(double mu, double link_power);
double linkInv(double eta, double link_power);
double mu_eta(double eta, double link_power) ;

void smat_multiply(double **A, double **B, double **C, int n) ;
void smat_inverse(double **A, double **B, int n) ;
