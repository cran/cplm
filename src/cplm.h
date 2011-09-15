//  
//  Header for functions used in cplm
//       Wayne Zhang
//

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <R_ext/Applic.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include "arms.h"


#define EPS (0.001)  // used in numerical gradient
#define SCALE (-1)   // used in optim

// struct to store data and parameters 
#ifndef DA_PARM_STRUCT 
#define DA_PARM_STRUCT

typedef struct {
  int nO ;              // number of Observations
  int nP ;              // # positive observations
  int nB ;              // # model cofficients
  int k ;               // index used internally by various functions
  int *ygt0 ;           // row index of positive values 
  double **X ;          // design matrix
  double *Y ;           // reponse variable 
  double *offset ;      // vector of offsets
  double *weights ;     // prior weights
  double *beta ;        // model coefficients
  double phi ;          // dispersion parameter
  double p ;            // index parameter
  double link_power ;   // power of link function, as in tweedie
  double lambda ;	// original mean for the truncated Poisson proposal	
}  da_parm;

#endif

// memory allocation utilities
double * dvect(int n) ;
double ** dmatrix(int nr, int nc) ;
int * ivect(int n) ;
int ** imatrix(int nr, int nc) ;

// simple computation
double dcumsum(double *x, int n) ;
int icumsum(int *x, int n) ;
double dcumwsum(double *x, double *w, int n) ;
double icumwsum(int *x, double *w, int n) ;
double norm (double *x, int n);
double dist (double *x, double *y, int n);
void smat_multiply(double **A, double **B, double **C, int n) ;
void smat_inverse(double **A, double **B, int n) ;

// univariate optimization 
void lbfgsbU(double *x, double lower, double upper, double *val,
             optimfn fn, optimgr gr, void *ex,  int *conv) ;

// glm related computation
double varFun (double mu, double p);
double linkFun(double mu, double link_power);
double linkInv(double eta, double link_power);
double mu_eta(double eta, double link_power) ;

// R list 
SEXP getListElement (SEXP list, char *str);

// cplm
double cplm_post_latT(double x, double y, double phi, double p) ; 
double cplm_lambda_tpois(double y, double phi, double p) ;
// function to simulate from 0 truncated poisson
int cplm_rtpois(double lambda ) ;
// function to compute the density of 0-truncated poisson on log scale
double cplm_dtpois(double x, double lambda) ;
void cplm_rlatT_reject (int nS, int *ans, da_parm *dap) ;

void cplm_tw_glm(double **x, double *y, double *off, double *wts, 
		 double *beta, double vp, double lp, 
		 double eps, int n, int p, int niter,
		 double  *xw, double *beta_old, double *zw,
		 double *qraux, double *work, double *qty, int *pivot) ;
void cpglm_fitted(double *beta, da_parm *da, double *mu,
		  double *eta) ;
double cplm_llikS(double *mu, double phi, double p,
		      int *simT, da_parm *dap);
