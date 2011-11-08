/************************************************************/
/*     Utility functions common to cplm programs            */
/*              Author:  Wayne Zhang                        */
/*            actuary_zhang@hotmail.com                     */
/************************************************************/


/**
 * @file utilities.c
 * @brief Utility functions common to cplm programs
 * @author Wayne Zhang                         
*/

#include "cplm.h"


/************************************************************/
/*          Memory allocation utility functions             */
/************************************************************/

/**
 * allocate memory for a double vector of size n
 *
 * @param n length of the desire vector
 *
 * @return a double pointer
 */
double * dvect(int n){
    return (double *)R_alloc(n, sizeof(double));
}


/**
 * allocate memory for an int vector of size n
 *
 * @param n length of the desire vector
 *
 * @return a int pointer
 */
int * ivect(int n){
    return (int *)R_alloc(n, sizeof(int));
}

/**
 * allocate memory for a double matrix of nr times nc
 *
 * @param nr number of rows
 * @param nc number of columns
 *
 * @return a 2d array
 */
double ** dmatrix(int nr, int nc)
{
    int i;
    double **m;
    m = (double **) R_alloc(nr, sizeof(double *));
    for (i = 0; i < nr; i++)
	m[i] = (double*) R_alloc(nc, sizeof(double));
    return m;
}

/**
 * allocate memory for an integer matrix of nr times nc
 *
 * @param nr number of rows
 * @param nc number of columns
 *
 * @return a 2d array
 */
int ** imatrix(int nr, int nc)
{
    int i;
    int **m;
    m = (int **) R_alloc(nr, sizeof(int *));
    for (i = 0; i < nr; i++)
	m[i] = (int*) R_alloc(nc, sizeof(int));
    return m;
}

/************************************************************/
/*               Simple arithmetic utility                  */
/************************************************************/

/**
 * cumulative sum of a vector of double 
 *
 * @param x double vector to be summed
 * @param n number of elements
 *
 * @return cumulative sum
 */
double dcumsum(double *x, int n){
    int i;
    double sm=0 ;
    for (i=0;i<n;i++)
        sm +=x[i] ;
    return sm ;
}

/**
 * cumulative sum of a vector of integer 
 *
 * @param x  vector to be summed
 * @param n number of elements
 *
 * @return cumulative sum
 */
int icumsum(int *x, int n){
    int i;
    int sm=0 ;
    for (i=0;i<n;i++)
        sm +=x[i] ;
    return sm ;
}

/**
 * weighted cumulative sum of a double vector by weight w
 *
 * @param x double vector to be summed
 * @param w weights
 * @param n number of elements
 *
 * @return cumulative sum
 */
double dcumwsum(double *x, double *w, int n){
    int i;
    double sm=0 ;
    for (i=0;i<n;i++)
        sm +=x[i]*w[i] ;
    return sm ;
}


/**
 * weighted cumulative sum of an int vector by weight w
 *
 * @param x vector to be summed
 * @param w weights
 * @param n number of elements
 *
 * @return cumulative sum
 */
double icumwsum(int *x, double *w, int n){
    int i;
    double sm=0 ;
    for (i=0;i<n;i++)
        sm +=x[i]*w[i] ;
    return sm ;
}


/**
 * take the norm of a vector
 *
 * @param x a double vector
 * @param n length of the vector
 *
 * @return norm of the vector
 */
double norm (double *x, int n){
    int i ;
    double nm=0 ;
    for (i=0;i<n;i++)
        nm += x[i] * x[i] ;
    return sqrt(nm) ;
}


/**
 * norm distance between two vectors
 *
 * @param x a double vector
 * @param y a double vector
 * @param n length of the vector
 *
 * @return norm distance of the two vectors
 */
double dist (double *x, double *y, int n){
    int i ;
    double nm=0 ;
    for (i=0;i<n;i++)
        nm += (x[i] - y[i])*(x[i] - y[i]) ;
    return sqrt(nm) ;
}

/**
 * get the max value of a double vector
 *
 * @param x a double vector
 * @param n length of the vector
 *
 * @return max value
 */
double dmax (double *x, int n){
    double ans = x[0] ;
    if (n>1){
        for (int i=1; i<n; i++)
            if (x[i]>ans) ans = x[i] ; 
    }
    return ans ;
}

/**
 * get the max value of an int vector
 *
 * @param x an int vector
 * @param n length of the vector
 *
 * @return max value
 */
int imax (int *x, int n){
    int ans = x[0] ;
    if (n>1){
        for (int i=1; i<n; i++)
            if (x[i]>ans) ans = x[i] ; 
    }
    return ans ;
}


/************************************************************/
/*                 Matrix computations                      */
/************************************************************/


/**
 * Multiply a matrix and a vector 
 *
 * @param trans transpose of matrix?
 * @param m row count of matrix
 * @param n column count of matrix
 * @param A input matrix
 * @param x input vector
 * @param out output vector 
 *
 */
void mult_mv(char *trans, int m, int n, double *A,
             double *x, double *out){
    double one = 1, zero = 0 ;
    int incx = 1;
    F77_CALL(dgemv)(trans, &m, &n, &one, A, &m, x, &incx,
                    &zero, out, &incx) ;
}


/**
 * compute t(x) * x
 *
 * @param m row dimension of the matrix
 * @param n column dimension of the matrix
 * @param x the input matrix  
 * @param out output results
 *
 */

void mult_xtx(int m, int n, double *x, double *out){
    double alpha=1.0, beta=0,
        *x2 = Calloc(m*n,double);
    Memcpy(x2, x, m*n) ;
    F77_CALL(dgemm)("T", "N", &n, &n, &m, &alpha, x2, &m,
                    x, &m, &beta, out, &n) ;
    Free(x2) ;
}

/**
 * compute the lower cholesky factor
 *
 * @param d dimension of the matrix
 * @param v input matrix
 * @param iv output cholesky factor
 *
 */
void chol(int d, double *v, double *iv){    
    int info;
    // cholesky factor of v
    Memcpy(iv, v, d*d) ;   
    F77_CALL(dpotrf)("L",&d,iv,&d,&info) ;
    if (info!=0) 
        error(_("Error %d in Cholesky decomposition."), info) ;   
}

/**
 * invert a positive symmetric matrix 
 *
 * @param d dimension of the matrix
 * @param v input matrix
 * @param iv output inverse of the matrix
 *
 */
void solve_po(int d, double *v, double *iv){    
    int info, i, j;
    // cholesky factor of v
    chol(d, v, iv) ;
    // compute inverse    
    F77_CALL(dpotri)("L",&d,iv,&d,&info) ;    
    if (info!=0) 
        error(_("Error %d in inverting matrix."), info) ;
    // fill upper triangle 
    for (i=0;i<d-1;i++){
        for (j=i+1;j<d;j++)
            iv[j*d+i] = iv[i*d+j] ;
    }    
}



/************************************************************/
/*               Optimization utility function              */
/************************************************************/

/**
 * optimation using the lbfgsb algorithm. This is a modification
 * of R's function lbfgsb, where memory allocation is changed
 * to Calloc and Free
 */
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
	error(_("REPORT must be > 0 (method = \"L-BFGS-B\")"));
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
		error(_("L-BFGS-B needs finite values of 'fn'"));
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
 *  wrapper for the univariate lbfgsb function:
 *
 * @param x the value of the parameter
 * @param lower the lower bound
 * @param upper the upper bound
 * @param val the value of the function
 * @param fn function to be minimized
 * @param gr derivative of function
 * @param ex struct to be passed to fn and gr
 * @param conv converged or not?
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

/**
 * extract element from a list
 *
 * @param list a R list
 * @param str the name of the element to be extracted
 *
 * @return a SEXP pointer to the required list element
 *
 **/
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



/************************************************************/
/*              cplm intermediate stats update              */
/************************************************************/

/**
 * compute variance function 
 *
 * @param mu mean 
 * @param p  index parameter
 *
 * @return variance function
 */ 
double varFun (double mu, double p){
    return pow(mu,p) ;
}


/**
 * get linear predictors 
 *
 * @param mu mean vector
 * @param link_power  link power
 *
 * @return linear predictors
 */ 
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


/**
 * get mean from linear predictors 
 *
 * @param eta linear predictors 
 * @param link_power  link power
 *
 * @return mean
 */ 
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

/**
 * derivative d mu / d eta
 *
 * @param eta linear predictors 
 * @param link_power  link power
 *
 * @return d(mu)/d(eta)
 */ 
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

/**
 * Compute the linear predictors given the design matrix and beta
 *
 * @param eta pointer to the vector of linear predictors
 * @param nO number of obs
 * @param nB number of coefficients
 * @param X pointer to the design matrix (in long vector representation)
 * @param beta pointer to the vector of coefficients
 * @param offset pointer to the vector of offsets
 */   

void cpglm_eta(double *eta, int nO, int nB, double *X, 
              double *beta,  double *offset){    
    // eta = X %*% beta
    mult_mv("N", nO, nB, X, beta, eta) ;
    // eta = eta + offset
    for (int i=0;i<nO;i++)
        eta[i] += offset[i] ;
}

/**
 * update the expected value and d(mu)/d(eta)
 *
 * @param mu pointer to the vector of expected values
 * @param muEta pointer to the vector of d(mu)/d(eta)
 * @param eta pointer to the vector of linear predictors
 * @param nO number of observations
 * @param link_power the power of the link function
 *
 */
void cplm_mu_eta(double *mu, double* muEta, int nO, 
                   double* eta, double link_power){
  int i ;
  for (i=0;i<nO;i++){
    mu[i] = linkInv(eta[i], link_power);
    if (muEta != NULL)
        muEta[i] = mu_eta(eta[i], link_power) ;
  }
}

/**
 * update the eta, mu and d(mu)/d(eta) in cpglm
 *
 * @param eta pointer to the vector of linear predictors
 * @param mu pointer to the vector of expected values
 * @param muEta pointer to the vector of d(mu)/d(eta)
 * @param dap pointer to a da_parm struct
 *
 */

void cpglm_fitted(double *eta, double *mu, double *muEta,
                 da_parm *dap){
    cpglm_eta(eta, dap->nO, dap->nB, dap->X, dap->beta,
             dap->offset) ;
    cplm_mu_eta(mu, muEta, dap->nO, eta, dap->link_power);
}


/**
 * update the power variance function
 *
 * @param var pointer to the vector of variances
 * @param mu pointer to the vector of mean values
 * @param n number of observations
 * @param p the index parameter
 *
 */
void cplm_varFun(double* var, double* mu, int n, double p){
  int i ;
  for (i=0;i<n;i++)
    var[i] = varFun(mu[i], p) ;
}


/************************************************************/
/*                Compute sample statistics                 */
/************************************************************/

/**
 * Compute sample mean
 *
 * @param n number of samples
 * @param x samples in long vector 
 *
 * @return mean 
 */
double mean(int n, double *x){
    return dcumsum(x, n)/n ;
}

/**
 * Compute sample variance for a univariate variable
 *
 * @param n number of samples
 * @param x samples in long vector 
 * @param ans pointer to store computed variance
 *
 */
void cov1(int n, double *x, double *ans){
    int i;
    double m = mean(n, x) ;
    ans[0] = 0;
    for (i=0;i<n;i++)
        ans[0] += (x[i]-m)*(x[i]-m);
    ans[0] /= n - 1.0  ;
}

/**
 * Compute sample covariance matrix 
 *
 * @param n number of samples
 * @param p number of variables (columns), p >2
 * @param x samples in long vector 
 * @param ans vector to store computed covariance matrix
 *
 */
void cov2(int n, int p, double *x, double *ans){
    int i;
    double *one = Calloc(n*n, double),
        *x2 = Calloc(n*p, double),
        *x3 = Calloc(n*p, double);
    double alpha = -1.0/n, beta = 1.0, beta2=0;

    // subtract mean    
    for (i=0;i<n*n;i++)
        one[i] = 1.0 ;
    Memcpy(x2, x, n*p) ;
    Memcpy(x3, x, n*p); 
    F77_CALL(dgemm)("N","N",&n,&p,&n, &alpha, one,
                    &n, x2, &n, &beta, x3, &n);
    Memcpy(x2, x3, n*p) ;
    AZERO(ans,p*p) ;
    
    // compute covariance 
    F77_CALL(dgemm)("T","N",&p,&p,&n, &beta, x2,
                    &n, x3, &n, &beta2, ans, &p);
    for (i=0;i<p*p;i++)
        ans[i] /= n-1 ;
    Free(one) ;
    Free(x2) ;
    Free(x3);
}

/**
 * Compute sample covariance matrix 
 *
 * @param n number of samples
 * @param p number of variables (columns)
 * @param x samples in long vector 
 * @param ans vector to store computed covariance matrix
 *
 */
void cov(int n, int p, double *x, double *ans){
    if (p==1)
        cov1(n, x, ans) ;
    else
        cov2(n, p, x, ans) ;
}



/************************************************************/
/*     Distribution and simulation related utilities        */
/************************************************************/

/**
 * simulation of multivariate normal
 *
 * @param d dimension
 * @param m mean vector
 * @param v positive-definite covarince matrix
 * @param s vector to store the simulated values
 *
 */

void rmvnorm(int d, double *m, double *v, double *s){
    int i, incx=1;
    double *lv = Calloc(d*d, double) ;
    GetRNGstate() ;
    // simulate d univariate normal r.v.s
    for (i=0;i<d;i++)
        s[i] = rnorm(0,1) ;    
    PutRNGstate() ;    
    // cholesky factor of v
    chol(d, v, lv) ;
    // scale and shift univariate normal r.v.s
    F77_CALL(dtrmv)("L","N","N",&d,lv,&d,s,&incx) ;
    for (i=0;i<d;i++)
        s[i] += m[i] ;    
    Free(lv) ;    
}

/**
 * return the exponent of a multivariate normal
 *
 * @param d dimension of the matrix
 * @param x sample vector 
 * @param m mean vector
 * @param iv inverse of the covariance matrix 
 *
 * @return exponent of MVN
 */
double dmvnorm(int d, double *x, double *m, double *iv){
    int i ;
    double ep=0, *dx = Alloca(d,double), *tmp = Alloca(d,double) ;
    R_CheckStack() ;
    for (i=0;i<d;i++)
        dx[i] = (m==NULL) ? x[i] : x[i] - m[i];
    mult_mv("N",d,d,iv,dx,tmp) ;
    for (i=0;i<d;i++)
        ep += dx[i]*tmp[i] ;
    ep *= -0.5 ;
    return ep ;
}

/**
 * random walk metropolis sampling for a vector of parameters 
 * of length d using multivariate normal proposal
 *
 * @param d dimension of the parameter
 * @param m current values of the parameter (also mean vector
 *      in the multivariate Normal)
 * @param v covariance matrix in the proposal distribution
 * @param sn simulated new vector 
 * @param myfunc user specified function to compute log posterior 
 * @param data struct used in myfunc
 *
 * @return  a 0-1 integer: 0 means not accept and 1 accept
 *
 */
int metrop_mvnorm_rw(int d, double *m, double *v, double *sn, 
		     double (*myfunc)(double *x, void *data), 
		     void *data){
    double A ;
    rmvnorm(d, m, v, sn) ;
    // determine if accept the sample
    A = exp(myfunc(sn,data)-myfunc(m,data) ) ;
    if (A<1 && runif(0,1)>=A){
        Memcpy(sn, m, d) ;
        return 0 ;
    }
    else return 1 ;
  
}


/**
 * Simulate truncated Normal r.v.s.
 *
 * @param m mean of the untrucated normal
 * @param sd standard deviation of the untrucated normal
 * @param lb left bound of the trucated normal
 * @param rb right bound of the trucated normal
 *
 * @return one simulated truncated normal 
 */

double cplm_rtnorm(double m, double sd, double lb, double rb){
    double u = runif(R_FINITE(lb)? pnorm(lb,m,sd,1,0): 0,
                     R_FINITE(rb)? pnorm(rb,m,sd,1,0): 1);
    return qnorm(u,m,sd,1,0) ;
}


/**
 * compute log density of truncated normal
 *
 * @param x the point at which to compute the log density
 * @param m mean of the untrucated normal
 * @param sd standard deviation of the untrucated normal
 * @param lb left bound of the trucated normal
 * @param rb right bound of the trucated normal
 *
 * @return log density at the point x
 */
double cplm_dtnorm(double x, double m, double sd, double lb, double rb){
    double c = R_FINITE(rb)? pnorm(rb,m,sd,1,0): 1 
        - R_FINITE(lb)? pnorm(lb,m,sd,1,0): 0 ;                     
    return dnorm(x,m,sd,1)- log(c) ;
}


/**
 * RW Metropolis update using trucated Normal 
 *
 * @param m mean of the untrucated normal
 * @param sd standard deviation of the untrucated normal
 * @param lb left bound of the trucated normal
 * @param rb right bound of the trucated normal
 * @param sn pointer to store simulated value
 * @param myfunc user specified function to compute log posterior 
 * @param data struct used in myfunc
 *
 * @return  a 0-1 integer: 0 means not accept and 1 accept
 */
int metrop_tnorm_rw( double m, double sd, double lb, double rb, double *sn, 
		     double (*myfunc)(double x, void *data), 
		     void *data){
    double A ;
    *sn = cplm_rtnorm(m, sd, lb, rb) ;
    // determine if accept the sample
    A = exp(myfunc(*sn,data)-myfunc(m,data)+
            cplm_dtnorm(m,*sn,sd,lb,rb)-
            cplm_dtnorm(*sn,m,sd,lb,rb)) ;
    if (A<1 && runif(0,1)>=A){ 
        *sn = m ;
        return 0 ;
    }
    else return 1 ;  
}

/**
 * compute log density for tweedie with positive response
 *
 * @param y  response
 * @param mu mean
 * @param phi scale parameter
 * @param p index parameter
 *
 * @return log density
 */

static double dtweedie2(double y, double mu,
                 double phi, double p){
    double a, a1, logz, drop = 37, jmax, j, cc, wmax, estlogw;
    double wm = -1.0E16, sum_ww = 0, *ww, ld;
    int k, lo_j, hi_j ;
    
    a= (2-p)/(1-p) ;
    a1 = 1 - a ;
    logz = -a *log(y) + a*log(p-1)- a1*log(phi)-log(2-p);
    jmax = pow(y,2-p)/(phi*(2-p)) ;

    jmax = fmax2(1.0, jmax) ;
    j = jmax ;
    cc = logz + a1 + a *log(-a) ;
    wmax = a1 * jmax ;
    estlogw = wmax ;
    while (estlogw > (wmax - drop)) {
        j += 2.0 ;
        estlogw = j * (cc - a1 * log(j)) ;
    }
    hi_j = ceil(j) ;
    j = jmax ;
    estlogw = wmax ;
    while ((estlogw > (wmax - drop)) && (j >= 2)) {
        j = fmax2(1, j - 2) ;
        estlogw = j * (cc - a1 * log(j)) ;
    }
    lo_j = imax2(1, floor(j)) ;
    ww = Calloc(hi_j-lo_j+1, double) ;
    for (k=lo_j;k<hi_j+1;k++){
        ww[k-lo_j] = k*logz - lgamma(1+k) - lgamma(-a*k) ;
        wm = fmax2(wm, ww[k-lo_j]) ;
    }
    for (k=lo_j;k<hi_j+1;k++)
        sum_ww += exp(ww[k-lo_j]-wm) ;

    ld = -y/(phi * (p - 1) * pow(mu,p - 1)) -
        (pow(mu, 2 - p)/(phi * (2 - p)))- log(y) +
        log(sum_ww) + wm  ;
    Free(ww) ;
    return ld ;
            
}


/**
 * compute log density for tweeide 
 *
 * @param y  response
 * @param mu mean
 * @param phi scale parameter
 * @param p index parameter
 *
 * @return log density
 */
double dtweedie(double y, double mu,
                 double phi, double p){
    return y ? dtweedie2(y,mu, phi,p):
        (-pow(mu,2 - p)/(phi * (2 - p))) ;
         
}

/**
 * compute -2 logliklihood of tweedie
 *
 * @param n number of obs
 * @param y vector of response
 * @param mu vector of means
 * @param phi scale parameter
 * @param p index parameter
 *
 * @return -2 loglik
 */

double dl2tweedie(int n, double *y, double *mu,
                  double phi, double p){
    int i ;
    double ans = 0;
    for (i=0;i<n;i++)
        ans += dtweedie(y[i],mu[i],phi,p) ;
    ans *= -2 ;
    return ans ;

}

/**
 * Simulate the Cholesky factor of a standardized Wishart variate with
 * dimension p and nu degrees of freedom.
 *
 * @param nu degrees of freedom
 * @param p dimension of the Wishart distribution
 * @param upper if 0 the result is lower triangular, otherwise upper
                triangular
 * @param ans array of size p * p to hold the result
 *
 * @return ans
 */
static double *std_rWishart_factor(double nu, int p, int upper, double ans[])
{
    int pp1 = p + 1;

    if (nu < (double) p || p <= 0)
	error("inconsistent degrees of freedom and dimension");

    AZERO(ans, p * p);
    for (int j = 0; j < p; j++) {	/* jth column */
	ans[j * pp1] = sqrt(rchisq(nu - (double) j));
	for (int i = 0; i < j; i++) {
	    int uind = i + j * p, /* upper triangle index */
		lind = j + i * p; /* lower triangle index */
	    ans[(upper ? uind : lind)] = norm_rand();
	    ans[(upper ? lind : uind)] = 0;
	}
    }
    return ans;
}

/**
 * Simulate a sample of random matrix from a Wishart distribution
 *
 * @param d row (=column) dimension of the matrix 
 * @param nu Degrees of freedom
 * @param scal Positive-definite scale matrix
 * @param out simulated matrix (d*d)
 *
 */
void rwishart(int d, double nu, double *scal, double *out)
{
    int  info,  psqr;
    double *scCp, *tmp, one = 1, zero = 0;

    psqr = d*d;
    tmp = Calloc(psqr, double);
    scCp = Calloc(psqr, double);

    Memcpy(scCp, scal, psqr);
    AZERO(tmp, psqr);
    F77_CALL(dpotrf)("U", &d, scCp, &d, &info);
    if (info)
	error(_("scal matrix is not positive-definite"));
    GetRNGstate();    
    std_rWishart_factor(nu, d, 1, tmp);
    F77_CALL(dtrmm)("R", "U", "N", "N", &d, &d,
			&one, scCp, &d, tmp, &d);
    F77_CALL(dsyrk)("U", "T", &d, &d, &one, tmp, &d,
			&zero, out, &d);
    for (int i = 1; i < d; i++){
        for (int k = 0; k < i; k++)
            out[i + k * d] = out[k + i * d];
    }
    PutRNGstate();
    Free(tmp) ;
    Free(scCp) ;
}
