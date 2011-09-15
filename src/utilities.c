//  
//  Some utility functions used in cplm
//       Wayne Zhang
//

#include "cplm.h"

// allocate memory for a double vector of size n
double * dvect(int n)
{
    return (double *)R_alloc(n, sizeof(double));
}

// allocate memory for a double vector of size n
int * ivect(int n)
{
    return (int *)R_alloc(n, sizeof(int));
}

// allocate memory for a double matrix of nr \times nc
double ** dmatrix(int nr, int nc)
{
    int i;
    double **m;
    m = (double **) R_alloc(nr, sizeof(double *));
    for (i = 0; i < nr; i++)
	m[i] = (double*) R_alloc(nc, sizeof(double));
    return m;
}

// allocate memory for an integer matrix of nr \times nc
int ** imatrix(int nr, int nc)
{
    int i;
    int **m;
    m = (int **) R_alloc(nr, sizeof(int *));
    for (i = 0; i < nr; i++)
	m[i] = (int*) R_alloc(nc, sizeof(int));
    return m;
}

// cumulative sum of a vector of double 
double dcumsum(double *x, int n){
    int i;
    double sm=0 ;
    for (i=0;i<n;i++)
        sm +=x[i] ;
    return sm ;
}

// cumulative sum of a vector of integer 
int icumsum(int *x, int n){
    int i;
    int sm=0 ;
    for (i=0;i<n;i++)
        sm +=x[i] ;
    return sm ;
}

// weight cumulative sum of a double vector by weight w
double dcumwsum(double *x, double *w, int n){
    int i;
    double sm=0 ;
    for (i=0;i<n;i++)
        sm +=x[i]*w[i] ;
    return sm ;
}

// weight cumulative sum of an int vector by weight w
double icumwsum(int *x, double *w, int n){
    int i;
    double sm=0 ;
    for (i=0;i<n;i++)
        sm +=x[i]*w[i] ;
    return sm ;
}


void lbfgsb2(int n, int m, double *x, double *l, double *u, int *nbd,
	    double *Fmin, optimfn fminfn, optimgr fmingr, int *fail,
	    void *ex, double factr, double pgtol,
	    int *fncount, int *grcount, int maxit, char *msg,
	    int trace, int nREPORT)
{
    char task[60];
    double f, *g, dsave[29], *wa;
    int tr = -1, iter = 0, *iwa, isave[44], lsave[4];

    if(n == 0) { /* not handled in setulb */
	*fncount = 1;
	*grcount = 0;
	*Fmin = fminfn(n, u, ex);
	strcpy(msg, "NOTHING TO DO");
	*fail = 0;
	return;
    }
    if (nREPORT <= 0)
	error(("REPORT must be > 0 (method = \"L-BFGS-B\")"));
    switch(trace) {
    case 2: tr = 0; break;
    case 3: tr = nREPORT; break;
    case 4: tr = 99; break;
    case 5: tr = 100; break;
    case 6: tr = 101; break;
    default: tr = -1; break;
    }

    *fail = 0;
    g = Calloc(n, double);
    /* this needs to be zeroed for snd in mainlb to be zeroed */
    wa = Calloc(2*m*n+4*n+11*m*m+8*m, double);
    iwa = Calloc(3*n, int);
    strcpy(task, "START");
    while(1) {
	/* Main workhorse setulb() from ../appl/lbfgsb.c : */
	setulb(n, m, x, l, u, nbd, &f, g, factr, &pgtol, wa, iwa, task,
	       tr, lsave, isave, dsave);
/*	Rprintf("in lbfgsb - %s\n", task);*/
	if (strncmp(task, "FG", 2) == 0) {
	    f = fminfn(n, x, ex);
	    if (!R_FINITE(f))
		error(("L-BFGS-B needs finite values of 'fn'"));
	    fmingr(n, x, g, ex);
	} else if (strncmp(task, "NEW_X", 5) == 0) {
	    if(trace == 1 && (iter % nREPORT == 0)) {
		Rprintf("iter %4d value %f\n", iter, f);
	    }
	    if (++iter > maxit) {
		*fail = 1;
		break;
	    }
	} else if (strncmp(task, "WARN", 4) == 0) {
	    *fail = 51;
	    break;
	} else if (strncmp(task, "CONV", 4) == 0) {
	    break;
	} else if (strncmp(task, "ERROR", 5) == 0) {
	    *fail = 52;
	    break;
	} else { /* some other condition that is not supposed to happen */
	    *fail = 52;
	    break;
	}
    }
    *Fmin = f;
    *fncount = *grcount = isave[33];
    if (trace) {
	Rprintf("final  value %f \n", *Fmin);
	if (iter < maxit && *fail == 0) Rprintf("converged\n");
	else Rprintf("stopped after %i iterations\n", iter);
    }
    strcpy(msg, task);
    Free(g) ;
    Free(wa) ;
    Free(iwa) ;
}


/*
 wrapper for the univariate lbfgsb function:
 - x, the value of the parameter
 - lower, the lower bound
 - upper, the upper bound
 - val, the value of the function
 - fn, function to be minimized
 - gr, derivative of function
 - ex, struct to be passed to fn and gr
*/

void lbfgsbU(double *x, double lower, double upper, double *val,
             optimfn fn, optimgr gr, void *ex,  int *conv){
    // nbd needed for lbfgsb 
    int nbd ;
    if (!R_FINITE(lower)) {
        if (!R_FINITE(upper)) nbd = 0; else nbd = 3;
    } else {
        if (!R_FINITE(upper)) nbd = 1; else nbd = 2;
    }

    // default the parameters needed for lbfgsb 
    double factr=1E7, pgtol=0 ;
    int fncount, grcount;
    int lmm=5, maxit=1000, trace=0, nREPORT=10 ;
    char msg[60] ;

    lbfgsb2(1,lmm, x, &lower, &upper, &nbd, val, fn, gr,
           conv, (void *)ex, factr, pgtol, &fncount, &grcount,
           maxit,msg ,trace, nREPORT) ;
} 


// extract element from a list
SEXP getListElement (SEXP list, char *str)
{
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;
  for (i = 0; i < length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  return elmt;
}



// take the norm of a vector
double norm (double *x, int n){
    int i ;
    double nm=0 ;
    for (i=0;i<n;i++)
        nm += x[i] * x[i] ;
    return sqrt(nm) ;
}

// norm distance between two vectors
double dist (double *x, double *y, int n){
    int i ;
    double nm=0 ;
    for (i=0;i<n;i++)
        nm += (x[i] - y[i])*(x[i] - y[i]) ;
    return sqrt(nm) ;
}

// variance function 
double varFun (double mu, double p){
    if (p==0)
        return 0.0 ;
    if (p==1)
        return mu ;
    if (p==2)
        return mu *mu ;
    else 
        return pow(mu,p) ;
}

// link function 
double linkFun(double mu, double link_power){
    if (link_power==0) 
        return log(mu);
    else if (link_power==1)
        return mu ;
    else if (link_power==2) 
        return mu * mu ;
    else if (link_power==-1)
        return 1.0/mu ;
    else 
        return pow(mu,link_power) ;
}

// inverse link function
double linkInv(double eta, double link_power){
    if (link_power==0) 
        return exp(eta);
    else if (link_power==1)
        return eta ;
    else if (link_power==2) 
        return sqrt(eta) ;
    else if (link_power==-1)
        return 1.0/eta ;
    else 
        return pow(eta,1/link_power) ;
}

// derivative d mu / d eta
double mu_eta(double eta, double link_power){
    if (link_power==0) 
        return exp(eta);
    else if (link_power==1)
        return 1.0 ;
    else if (link_power==2) 
        return 1.0/(sqrt(eta)*2) ;
    else if (link_power==-1)
        return -1.0/(eta*eta) ;
    else 
        return pow(eta,1/link_power-1)/link_power ;
}

/*
// deviance residuals based on quasi-likelihood
double dev_resids(double y, double mu, double wts, double vp ){
    // ad hoc adjustment for y==0 
    double y1, theta, kappa, p1=1-vp, p2=2-vp ;
    y1 = (y==0) ? 1 : y;
    theta = (vp==1) ? (log(y1) -log(mu)) : ((pow(y,p1) - pow(mu,p1))/p1) ;
    kappa = (vp==2) ? (log(y1) - log(mu)) : ((pow(y,p2) - pow(mu,p2))/p2) ;
    return 2 * wts * (y*theta-kappa) ;
}
*/

// function to compute inverse of a general matrix
void smat_inverse(double **A, double **B, int n){
  int i, j, info, *ipiv, lwork=n*n;
  double *lA = Calloc(n*n, double), 
    *work = Calloc(n*n, double);
  ipiv = Calloc(n, int) ;
  for (i=0;i<n;i++){
    for (j=0;j<n;j++)
      lA[i+j*n] = A[i][j] ;
  }
  F77_CALL(dgetrf)(&n,&n,lA,&n,ipiv,&info) ;
  if (info!=0) 
    error("LU decomposition failed.") ;
  F77_CALL(dgetri)(&n,lA,&n,ipiv,work, &lwork, &info);
  if (info!=0) 
    error("Matrix inversion failed.") ;
  for (i=0;i<n;i++){
    for (j=0;j<n;j++)
       B[i][j] =lA[i+j*n];
  }
  Free(lA) ;
  Free(work) ;
  Free(ipiv) ;
}

// function to multiply two squared matrix 
void smat_multiply(double **A, double **B, double **C, int n){
  int i, j, k ;
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      C[i][j] = 0 ;
      for (k=0;k<n;k++)
	C[i][j] += A[i][k] * B[k][j] ;
    }
  }
}
