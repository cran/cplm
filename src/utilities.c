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
 * @param b pointer to the vector of random effects (unnormalized)
 * @param offset pointer to the vector of offsets
 */   

void cplm_eta(double *eta, int nO, int nB, double *X, 
              double *beta, double *b, double *offset){    
    double one= 1.0 ;
    int inc = 1, i;
    // initialize eta
    if (b==NULL) {
        AZERO(eta, nO) ;
    }
    else {
        Memcpy(eta, b, nO) ;
    }
    // eta = X %*% beta + b
    F77_CALL(dgemv)("N",&nO,&nB,&one,X,&nO,beta,&inc,&one,eta,&inc) ;
    // eta = eta + offset
    for (i=0;i<nO;i++)
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
    cplm_eta(eta, dap->nO, dap->nB, dap->X, dap->beta,
             (double *) NULL, dap->offset) ;
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



/**
 * simulation of multivariate normal
 *
 * @param d dimension
 * @param m mean vector
 * @param v positive-definite covarince matrix
 * @param s vector to store the simulated values
 *
 */
// FIXME: does it handle d==1?

void rmvnorm(int d, double *m, double *v, double *s){
    int i, info, incx=1;
    double *lv = Calloc(d*d, double) ;
    GetRNGstate() ;
    // simulate d univariate normal r.v.s
    for (i=0;i<d;i++)
        s[i] = rnorm(0,1) ;
    PutRNGstate() ;    
    // cholesky factor of v
    Memcpy(lv, v, d*d) ;
    F77_CALL(dpotrf)("L",&d,lv,&d,&info) ;
    if (info!=0) 
        error(_("Error %d in Cholesky decomposition."), info) ;        
    // scale and shift univariate normal r.v.s
    F77_CALL(dtrmv)("L","N","N",&d,lv,&d,s,&incx) ;
    for (i=0;i<d;i++)
        s[i] += m[i] ;    
    Free(lv) ;    
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
    rmvnorm(d,m, v, sn) ;
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
 * Compute sample covariance matrix 
 *
 * @param n number of samples
 * @param p number of variables (columns)
 * @param x samples in long vector 
 * @param ans vector to store computed covariance matrix
 *
 */

 void cplm_cov(int n, int p, double *x, double *ans){
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
    F77_CALL(dgemm)("N","N",&n,&p,&n, &alpha, one, &n, x2, &n, &beta, x3, &n);
    Memcpy(x2, x3, n*p) ;
    AZERO(ans,p*p) ;
    
    // compute covariance 
    F77_CALL(dgemm)("T","N",&p,&p,&n, &beta, x2,&n, x3, &n, &beta2, ans, &p);
    for (i=0;i<p*p;i++)
        ans[i] /= n-1 ;
    Free(one) ;
    Free(x2) ;
    Free(x3);
}

/**
 * compute log density for tweeide with positive response
 *
 * @param y  response
 * @param mu mean
 * @param phi scale parameter
 * @param p index parameter
 *
 * @return log density
 */

double dtweedie2(double y, double mu,
                 double phi, double p){
    double a, a1, logz, drop = 37, jmax, j, cc, wmax,
        estlogw, oldestlogw;
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
        oldestlogw = estlogw ;
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
