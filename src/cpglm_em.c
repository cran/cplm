/************************************************************/
/*   Function to implement the Monte Carlo EM algorithm     */
/*    in the Compound Poisson Generalized Linear Model      */
/*              Author:  Wayne Zhang                        */
/*            actuary_zhang@hotmail.com                     */
/*                   09/06/2011                             */
/************************************************************/

#include "utilities.h"

/* This is the struct needed for the MCEM glm 
   nO, nP, nS, nB: number of Obs, positive Obs, simulations, and parameters
   k: index used internally by various functions
   ygt0: row index of positive values 
   simT: 2d array of simulated latent Poisson variable
   X: design matrix
   Y: vector of all response  	       
   offset, weights: offset or weights in glm
   beta, phi, p: parameters in glm
   link_power: power of link function, as in tweedie
   lambda: original mean for the truncated Poisson proposal
   et: expected value of latent Poisson E(simT)
   iw, sw: importance weights and their sum
   mu, eta, mu1p, mueta: fitted expected value, linear predictor, 
        mu^(1-p) and d(mu)/d(eta) as used in computing gradient and hessian
*/

typedef struct {
  int nO, nP, nS, nB, 
    k, *ygt0, **simT; 
  double **X, *Y, *offset, *weights, 
    *beta, phi, p, link_power, 
    lambda, *et, *iw, sw, 
    *mu, *eta, *mu1p, *mueta;		
}  cpglm_str;


/* 
   Log posterior of T_i 
   - x, a point at which the density to be computed
   - y, the positive response value corresponding to the t_i
   - phi, scale parameter (if weights supplied, phi/weights)
   - p, index parameter
*/

double dLogPostT(double x, double y, double phi, double p){
  double p1=p-1, p2=2-p, ans;
  ans = x*(p2/p1 * (log(y) - log(p1)) - log(phi)/p1 -log(p2)) -
        lgamma(x*p2/p1) - lgamma(x+1) ;
  return ans ;
}

/*
  function to compute the mean of the proposal 0-truncated poisson 
  by finding the approximate mode
*/

double lambdaTruncPois(double y, double phi, double p){
  double lamdba ;
  lamdba = ceil(pow(y,2-p)/((2-p)*phi));
  return lamdba ;
}

// function to simulate from 0 truncated poisson
int rTruncPois(double lambda ){
  int sim=qpois(runif(dpois(0, lambda,0), 1), lambda, 1,0);
  return sim ;
}

// function to compute the density of 0-truncated poisson on log scale
double dLogTruncPois(int x, double lambda){
  double ld ;
  ld=x*log(lambda)- lambda -lgamma(x+1)-log(1-exp(-lambda)) ;
  return(ld) ;
}

/*
 compute difference between target and proposal (unnormalized)
 used to locate the point where max diff occurs
*/ 
static double sampleTDiff (double *par,  void *data){
    cpglm_str *da = data ;
    int k = da->k ;
    double phi2 = da->phi / da->weights[da->ygt0[k]] ;
    double ldiff = dLogPostT(par[0],da->Y[da->ygt0[k]],phi2, da->p)
        - dLogTruncPois(par[0],da->lambda) ;
    return ldiff ;
}

/*
  derivative of the sampleTDiff 
*/
static void sampleTGrad (double *par, double *gr,  void *data){
    cpglm_str *da = data ;
    int k= da->k ;
    double p1 = da->p -1 ,
        p2 = 2 - da->p ,
        phi2 = da->phi / da->weights[da->ygt0[k]] ;    
    *gr = p2/p1*log(da->Y[da->ygt0[k]]/p1)-log(phi2)/p1-log(p2)-
        log(da->lambda)-p2/p1*digamma(p2/p1*par[0]);
}

// the following two functions scale the above two by -1 for max
static double sampleTDiffOptim (int n, double *par,  void *data){    
  double ldiff = sampleTDiff(par, data) * SCALE ;
  return ldiff ;
}

static void sampleTGradOptim (int n, double *par, double *gr,  void *data){
  sampleTGrad(par,gr,data) ;
  *gr = *gr * SCALE ;
}

/*
  Rejections sampling of the latent variable T  *
  ans: return value is a length(ygt0) \times nS matrix of
     simulated values
*/
void sampleTRej (int **ans, void *data){
    cpglm_str *da = data ;    
    // dimensions of return matrix
    int nc = da->nS, nr = da->nP,
        i, j,  xtemp, accept, conv;
    double par=2.0, val, pb, u, xtemp2, phi2;
    GetRNGstate() ;
    
    for (i=0; i<nr;i++){
        da->k=i ;
        // adjust for prior weights
        phi2 = da->phi / da->weights[da->ygt0[i]] ;        
        da->lambda = lambdaTruncPois(da->Y[da->ygt0[i]],phi2, da->p) ;
        // find t that results in the largest difference
        // between target and proposal
        lbfgsbU(&par, 1, INFINITY, &val,sampleTDiffOptim, sampleTGradOptim,
                (void *)da,  &conv) ;
	val /= SCALE ;
        for (j=0;j<nc;j++){
            accept =0 ;
            while (accept==0) {
                xtemp =rTruncPois(da->lambda);
                // sampleTDiff only takes double!
                xtemp2 = xtemp ;
                pb = exp(sampleTDiff(&xtemp2,da)-val);
                u = runif(0,1);
                if (u<=pb){
                    ans[i][j] = xtemp ;
                    accept =1 ;
                }
            }
        }
    }
    PutRNGstate() ;
}



/* function to fit tweedie glm
  x, y: 2-d array of design matrix and response variable
  off, wts: offsets and weights
  beta: starting values of beta
  vp: variance power
  lp: link power
  eps: tolarence level
  n, p, niter: number of obs, parameters and max iterations
  xw: long version of x, multiplied by working weight
  zw: linearized y times working weight   
  beta_old: store interim beta
  qraux, work, qty, pivot: as required in dqrdc and dqrsl
 */
void tw_glm(double **x, double *y, double *off, double *wts, 
	 double *beta, double vp, double lp, 
	 double eps, int n, int p, int niter,
         double  *xw, double *beta_old, double *zw,
	 double *qraux, double *work, double *qty, int *pivot){
    
    int i,j, k, job, info, nf,
        nO=n, nB=p, nx=niter ;
    double  dum,  eta, mu, mu_eta_val,z, w;

    // iterative re-weighted least squares
    for (k=0;k<nx;k++){
        // compute eta
        for (i=0;i<nO;i++){    
            eta = 0 ;
            for (j=0;j<nB;j++)
                eta += x[i][j] * beta[j] ; 
            eta += off[i] ;
            mu = linkInv(eta,lp) ;
            mu_eta_val = mu_eta(eta,lp);
            // linearization to get new response z
            z = eta + (y[i] - mu) / mu_eta_val -off[i] ;
            // working weight in least square
            w = sqrt(wts[i]/varFun(mu,vp) ) * mu_eta_val;
            // xw and zw are the design matrix and response in weighted ls
            for (j=0;j<nB;j++)
                xw[i+j*nO] = x[i][j] * w ;
            zw[i] = z * w ;
        }
        // old beta 
        for (j=0;j<nB;j++)
            beta_old[j] = beta[j] ;
        // qr decomposition (stored in xw now)
        job = 0 ; // no pivoting
        F77_CALL(dqrdc)(xw, &nO, &nO, &nB, qraux, pivot, work, &job) ;
        // weighted least square fit
        job = 100; // only want to compute beta
                   // but qty is also computed automatically by dqrsl
        F77_CALL(dqrsl)(xw, &nO, &nO, &nB, qraux, zw,
                        &dum, qty, beta, &dum, &dum, &job, &info) ;

        // check if non-finite value is generated
        nf=0 ; 
        for (i=0;i<nB;i++)
            nf += R_FINITE(beta[i]) ;
        if (nf<nB)
            error("non-finite coefficients are found \n");
        // branch out if converged
        if (dist(beta, beta_old, nB)/ norm(beta,nB) < eps)
            break ;    
    }
    
    if (dist(beta, beta_old, nB)/ norm(beta,nB) >= eps)
        warning("Algorithm did not converge - max iteration numbers reached!\n");
    
}

/*
  function to compute mu and eta of cpglm
*/
void cpglm_fitted(double *beta, void *data, double *mu,
             double *eta){
  cpglm_str *da= data ;
  double xb ;
  int i, j ;

  for (i=0; i< da->nO; i++) {
        xb =0 ;
        for (j=0;j<da->nB;j++)
            xb += da->X[i][j] * beta[j] ;
	xb += da->offset[i] ;
        eta[i] = xb ;        
        mu[i] = linkInv(xb,da->link_power);
    }
}

/*
 compute joint likelihood for one Single realization 
 (k in data) of simT given mu, phi and p
*/
double dJLogLikS(double *mu, double phi, double p, void *data){
    cpglm_str *da = data ;
    int i, k ;
    double lp=0, p2=2-p, p1=p-1;

    for (i=0; i<da->nO; i++)
        lp += pow(mu[i],p2) * da->weights[i];
    lp /= (- phi*p2) ;
    for (i=0; i<da->nP; i++){
	k = da->ygt0[i] ;
        lp += - da->Y[k]*pow(mu[k],-p1)*da->weights[k] /(phi*p1) +
	  dLogPostT(da->simT[i][da->k],da->Y[k],phi/da->weights[k],p);
    }
    return lp ;
}

/*
 compute Expected joint likelihood 
*/
double dJLogLikE(double *mu, double phi, double p, void *data){
    cpglm_str *da = data ;
    int i, j, k;
    double lp=0, elgt, p1=p-1, p2=2-p;

    // taking expectation of dJLogLikS
    for (i=0; i<da->nO; i++)
        lp += pow(mu[i],p2) * da->weights[i];
    lp /= - phi*p2 ;
    for (i=0; i<da->nP; i++){
        k = da->ygt0[i] ;
        elgt =0 ; 
        for (j=0; j<da->nS;j++)
            elgt += lgammafn(da->simT[i][j]*p2/p1)* da->iw[j] ;         
        elgt /= da->sw ;
        lp += - da->Y[k]*pow(mu[k],-p1)*da->weights[i]/(phi*p1)- elgt
            + da->et[i]*(p2/p1*(log(da->Y[k])-log(p1))-log(phi/da->weights[k])/p1 - log(p2));
    }
    
    return lp ;
}

/*
  function to compute the (expected) log posterior of p
*/
double dLogPostP(double x, void *data){
    cpglm_str *da = data ;
    return dJLogLikE(da->mu, da->phi, x, da) ;
}

/*
  The following two functions are to compute the log posterior
  of p and its gradient, as required in lfgsbU
*/
static double dLogPostPOptim(int n, double *par, void *data){
    cpglm_str *da =data ;    
    double ans= dLogPostP(*par, da) * SCALE;
    return ans ;
}

// compute numerical derivative
static void dLogPostPGradOptim (int n, double *par, double *gr,  void *data){
    cpglm_str *da =data ;
    double par1 = *par +EPS, par2 =*par-EPS ;
    *gr=(dLogPostP(par1,da)-dLogPostP(par2,da))/(2*EPS) * SCALE;
}


/*
  Compute weights in importance sampling 
 */
void import_weight(double phi_new, double p_new, double phi_old,
		   double p_old, void *data, double *iw){
  cpglm_str *da = data ;
  int nS = da->nS, nP = da->nP ;
  int i, j ;
  double iw_max=-1.0E20 ;
   // compute weight
   for (j=0;j<nS;j++){
     iw[j] =0.0 ;
     for (i=0;i<nP;i++){
       iw[j] += dLogPostT(da->simT[i][j],da->Y[da->ygt0[i]], phi_new, p_new)
	 -dLogPostT(da->simT[i][j],da->Y[da->ygt0[i]], phi_old, p_old);
     }
     iw_max = fmax2(iw_max, iw[j]) ;    
   }
   for (j=0;j<nS;j++){
       iw[j] = exp(iw[j]-iw_max) ; //subtract the max to avoid numerical issue
   }
}

// initialize importance weight to unity
void unit_weight(int nS, double *iw, double *sw){
    int i ;
    for (i=0;i<nS;i++)
        iw[i] = 1.0 ;
    *sw = dcumsum(iw,nS) ;
}

// store parameter values to a matrix theta 
void store_theta(double *theta, cpglm_str *da){
    int i ;
    for (i=0;i< da->nB;i++)
        theta[i]= da->beta[i] ;
    theta[da->nB] = da->phi ;
    theta[da->nB+1] = da->p ;
}

/*
  compute derivative for Single L given one simulation of T:
   - beta and phi is computed analytically
   - p is computed numerically 
*/ 
void dJLogLikGradS (double *par, double *gr,  void *data){
    cpglm_str *da =data ;
    int i, j, k ;
    double phi=par[da->nB], p=par[da->nB+1], p2=2-p, p1=p-1;

    // update mu, eta, mu1p and mueta according to par
    cpglm_fitted(par, da, da->mu, da->eta);
    for (i=0;i<da->nO;i++){
      da->mu1p[i]= pow(da->mu[i],1-p) ;
      da->mueta[i] = mu_eta(da->eta[i],da->link_power) ;
    }

    // derivative for beta
    for (j=0; j<da->nB;j++){
        gr[j] = 0 ;
        for (i=0;i<da->nO;i++)
            gr[j] += - da->weights[i]* da->mu1p[i]
                * da->mueta[i]*da->X[i][j];
        for (i=0;i<da->nP;i++){
            k = da->ygt0[i] ;
            gr[j] += da->Y[k]* da->mu1p[k]* da->weights[k]
                /da->mu[k]* da->mueta[k]*da->X[k][j] ;
        }
        gr[j] /= phi ;
    }

    // derivative for phi
    gr[da->nB]=0 ;
    for (i=0;i<da->nO;i++)
        gr[da->nB] += da->weights[i]*da->mu1p[i]*da->mu[i];
    gr[da->nB] /= (phi*phi*p2);
    for (i=0;i<da->nP;i++) {
        k = da->ygt0[i] ;
        gr[da->nB] += da->weights[k]*da->Y[k]*da->mu1p[k]/(p1*phi*phi)
            - da->simT[i][da->k]/(p1*phi) ;
    }

    // derivative for p, using numerical solution
    double dp1, dp2 ;
    p += EPS ;
    dp1 = dJLogLikS(da->mu,phi,p,da) ;
    p -= 2*EPS ;
    dp2 = dJLogLikS(da->mu,phi,p,da) ;
    gr[da->nB+1]=(dp1-dp2)/(2*EPS);
}

/*
  compute  derivative for E(L), MC average across all T's
*/ 

void dJLogLikGradE (double *par, double *gr,  void *data){
    cpglm_str *da =data ;
    int i, j, k ;
    double phi=par[da->nB], p=par[da->nB+1];
    double p2=2-p, p1=p-1;

    // update mu, eta, mu1p and mueta according to theta
    cpglm_fitted(par, da, da->mu, da->eta);
    for (i=0;i<da->nO;i++){
      da->mu1p[i]= pow(da->mu[i],1-p) ;
      da->mueta[i] = mu_eta(da->eta[i],da->link_power) ;
    }
    
    // derivative for beta
    for (j=0; j<da->nB;j++){
        gr[j] = 0 ;
        for (i=0;i<da->nO;i++)
            gr[j] += - da->weights[i] * da->mu1p[i]
                * da->mueta[i] * da->X[i][j];
        for (i=0;i<da->nP;i++){
            k = da->ygt0[i] ;
            gr[j] += da->weights[k] * da->Y[k]* da->mu1p[k]/ da-> mu[k]
                * da->mueta[k] * da->X[k][j] ;
        }
        gr[j] *= (1/phi) ;
    }

    // derivative for phi
    gr[da->nB]=0 ;
    for (i=0;i<da->nO;i++)
        gr[da->nB] += da->weights[i]*da->mu1p[i]*da->mu[i];
    gr[da->nB] /= (phi*phi*p2) ;
    for (i=0;i<da->nP;i++) {
        k = da->ygt0[i] ;
        gr[da->nB] += da->weights[k]*da->Y[k]* da->mu1p[k]/(p1*phi*phi)
            - da->et[i]/(p1*phi) ;
    }

    // derivative for p, using numerical solution
    double dp1, dp2 ;
    p += EPS ;
    dp1 = dJLogLikE(da->mu, phi, p,da) ;
    p -= 2*EPS ;
    dp2 = dJLogLikE(da->mu, phi, p,da) ;
    gr[da->nB+1]=(dp1-dp2)/(2*EPS);
}

// compute hessian numerically
void dJLogLikHessE (double *par, double **hess,  void *data){
    cpglm_str *da =data ;
    int i, j ;    
    double *df1 = dvect(da->nB+2);
    double *df2 = dvect(da->nB+2);

    for (i = 0; i < da->nB+2; i++) {
	par[i] +=  EPS;
	dJLogLikGradE(par, df1, da);
	par[i] -=  2 * EPS;
        dJLogLikGradE(par, df2, da);
	for (j = 0; j < da->nB+2; j++)
	    hess[i][j] = (df1[j] - df2[j])/ (2*EPS) ;
	par[i] +=EPS;
    }
}

/*
  Function to compute hessian using Monte Carlo methods
  - par, all the parameters in model: beta, phi,p in order
  - data, cpglm_str
  - VL, variance of L
  - HQ, hessian of E(L)
  - VQ, approx variance of E(L)
*/

void hessGlmEst (double *par,  void *data,
                 double **VL, double **HQ, double **VQ){
    cpglm_str *da = data ;
    int nS=da->nS, nB = da->nB ;
    int i, j, k;
    double **ELLt = dmatrix(nB+2,nB+2) ; // E(L*L')
    double *EL = dvect(nB+2) ;
    double *grad = dvect(nB+2) ;

    //initialize EL and ELLt
    for (j=0;j<nB+2;j++){
        EL[j] = 0 ;
        for (k=0;k<nB+2;k++)
            ELLt[j][k] = 0 ;
    }

    // compute E(L) and E(L*L') 
    for (i=0;i<nS;i++){
        da->k=i ;
        dJLogLikGradS(par, grad, da);
        for (j=0;j<nB+2;j++){
            // compute E(L)
            EL[j] += grad[j]/nS ;
            for (k=0;k<nB+2;k++)
                // compute E(L*L')
                ELLt[j][k] += grad[j]*grad[k]/nS;
        }
    }

    // compute VL and VQ
    for (i=0;i<nB+2;i++){
      for (j=0;j<nB+2;j++){
            VL[i][j]=ELLt[i][j]-EL[i]*EL[j] ;
	    VQ[i][j]=ELLt[i][j]/nS ;	    
      }
    }

    // compute hessian of Q (HQ)
    dJLogLikHessE(par, HQ, da);            
}


/************************************************/
/*     Main function to fit compound Poisson    */
/*     GLM using Monte Carlo EM algorithm       */
/************************************************/
/*
X, Y: design matrix and response vector
ygt0: row index of positive values 
offset, weights: offset or weights in glm
beta, phi, p: parameters in glm
link_power: power of link function, as in tweedie
bound: a vector giving the lower and upper bound of p
init_size, sample_iter, max_iter, max_size:
  initial sample size, # iterations after which to use
  importance sampling, max # iterations and max sample size
epsilon1, epsilon2, alpha, ck: as used in stopping rule
fixed_size: whether sample size should be fixed
trace: if fitting info should be printed
beta_step: integer, number of iterations to skip for updating beta
*/

SEXP cpglm_em (SEXP X, SEXP Y, SEXP ygt0, SEXP offset, SEXP weights, 
		SEXP beta, SEXP phi, SEXP p, SEXP link_power, SEXP bound, 
               SEXP init_size, SEXP sample_iter, SEXP max_iter,
	       SEXP epsilon1, SEXP epsilon2, SEXP alpha, SEXP ck, 
	       SEXP fixed_size, SEXP trace, SEXP max_size, SEXP beta_step){

  // dimensions 
    int nO=LENGTH(offset),
        nP = LENGTH(ygt0),
        nS=INTEGER(init_size)[0],
        nB=LENGTH(beta),
        iter_glm=100; //max glm iteration numbers
    int i, j, k, ii, pivot, conv, iter=0;    
    cpglm_str *da ;
    // glm related variables
    double *xw, *beta_old, *zw, *qraux, *work, *qty, eps_glm =1.0E-8;
    // em iteration variables 
    double **theta, *diff, p1, p2, dpval, pnew, T=0, 
        C = qchisq(1-REAL(alpha)[0], 2, 1,0);
    // var-cov related variables
    double  **VQ, **HQ, **VL, **mc_var_inv, **HQ2, **VQ2, **mc_var_inv2;

    // allocate memory for data struct
    da = (cpglm_str *) R_alloc(1,sizeof(cpglm_str)) ;
    da->X = dmatrix(nO,nB) ;
    da->Y = dvect(nO) ;
    da->ygt0 = ivect(nP) ;
    da->offset = dvect(nO) ;
    da->weights = dvect(nO) ;
    da->beta= dvect(nB) ;
    da->mu = dvect(nO) ;
    da->eta = dvect(nO) ;
    da->mu1p = dvect(nO) ;
    da->mueta = dvect(nO) ;
    da->et = dvect(nP) ;    
    if (INTEGER(fixed_size)[0]==1){ 
      da->simT = imatrix(nP, nS) ;
      da->iw = dvect(nS) ;
    }
    else {// set column size to max possible to avoid resizing
      da->simT = imatrix(nP, INTEGER(max_size)[0]) ; 
      da->iw = dvect(INTEGER(max_size)[0]) ;
    }
    // initialize to unit weight
    unit_weight(nS, da->iw, &(da->sw));

    // for em iterations
    diff = dvect(nB+2) ;
    VL = dmatrix(nB+2, nB+2) ;
    HQ = dmatrix(nB+2, nB+2) ;
    VQ = dmatrix(nB+2, nB+2) ;
    mc_var_inv = dmatrix(2, 2) ;
    HQ2 =dmatrix(2,2) ;
    VQ2 =dmatrix(2,2) ;
    mc_var_inv2 = dmatrix(2,2) ;
    // row dimension as max possible iterations
    theta = dmatrix(INTEGER(max_size)[0],nB+2) ;
    
    //This block is used for fitting glm    
    xw = dvect(nO*nB) ;
    beta_old = dvect(nB) ;
    zw = dvect(nO);
    qraux = dvect(nB) ;
    work = dvect(nB) ;
    qty = dvect(nO) ;
    
    // fill in struct da
    da->nO = nO ;
    da->nP = nP;    
    da->nS = nS ;
    da->nB = nB ;
    for (i=0;i<nP;i++)
        // ygt0 from R does not start from 0
        da->ygt0[i]=INTEGER(ygt0)[i]-1 ;
    for (i=0;i<nO;i++){
        da->Y[i]=REAL(Y)[i] ;
        da->offset[i]=REAL(offset)[i] ;
        da->weights[i]=REAL(weights)[i] ;
        for (j=0;j<nB;j++)
            da->X[i][j] = REAL(X)[i+j*nO];
    }
    for (i=0;i<nB;i++)
        da->beta[i] = REAL(beta)[i] ;
    da->phi = REAL(phi)[0] ;
    da->p = REAL(p)[0];
    da->link_power = REAL(link_power)[0] ;
    cpglm_fitted(da->beta, da, da->mu, da->eta);
    // store theta to begin with
    iter =0 ;
    store_theta(theta[iter], da) ;
    
    for (ii=0;ii<INTEGER(max_iter)[0];ii++){
        R_CheckUserInterrupt() ;
        iter++ ;
        /********************************************/
        /*************   E-step	   ******************/
        /********************************************/
        
	/** Rejection sampling of T **/
	if (iter <=INTEGER(sample_iter)[0]){	  
	  sampleTRej(da->simT,da) ;       
	  if (INTEGER(fixed_size)[0]==0)
              unit_weight(da->nS, da->iw, &(da->sw));                            
	}
	/** Importance sampling of T **/
	else {  
	  import_weight(da->phi, da->p, 
			theta[INTEGER(sample_iter)[0]-1][nB],
			theta[INTEGER(sample_iter)[0]-1][nB+1],
			da, da->iw) ;
	  da->sw = dcumsum(da->iw, da->nS) ;
	}	  	  
	// compute weighted E(T_i)
	for (i=0;i<nP;i++)
	  da->et[i]=icumwsum(da->simT[i],da->iw,da->nS)/da->sw;
        
  	/********************************************/
        /****************   M-step    ***************/	
        /*********************************************/
        
        /** update p **/        
	pnew = da->p ;
	lbfgsbU(&pnew, REAL(bound)[0], REAL(bound)[1], &dpval,
		dLogPostPOptim, dLogPostPGradOptim, (void *)da,  &conv) ;
	// if converged, set to new value; o/w, retain
	if (conv == 0) 
	  da->p = pnew ;

        R_CheckUserInterrupt() ;
	/** GLM update of beta **/
        // do this every (beta_step) iterations as beta varies little
        if ((iter-1)%INTEGER(beta_step)[0] ==0){
            tw_glm(da->X, da->Y, da->offset, da->weights, da->beta,
                   da->p, da->link_power, eps_glm, nO, nB, iter_glm,
                   xw, beta_old, zw, qraux, work, qty, &pivot);
            // compute fitted values and related quantities 
            cpglm_fitted(da->beta, da, da->mu, da->eta);
            for (i=0;i<nO;i++)
                da->mueta[i] = mu_eta(da->eta[i],da->link_power) ;            
        }
        for (i=0;i<nO;i++)
            da->mu1p[i]= pow(da->mu[i],1-da->p) ;

        /** update phi **/
	p1=da->p-1, p2=2-da->p;
        da->phi = 0;
	for (i=0;i<nO;i++)
	  da->phi += pow(da->mu[i],p2)*da->weights[i] ;
	da->phi *= p1/p2 ;
	for (i=0;i<nP;i++){
	  k = da->ygt0[i] ;
	  da->phi += da->Y[k]*da->weights[k]*pow(da->mu[k],-p1);
	}
	da->phi /= dcumsum(da->et,nP) ;
        
        // store results 
        store_theta(theta[iter], da) ;
                
        // print out iteration info if necessary  
        if (INTEGER(trace)[0]==1){
            if (iter==INTEGER(sample_iter)[0]+1)
                Rprintf("Importance sampling begins...\n") ;
            Rprintf("Iteration: %d, phi: %6.3f, p: %6.4f",iter,da->phi, da->p) ;        
        }
        // determine if stopping rule is reached 
        for (i=0;i<nB+2;i++)
	  diff[i] = fabs(theta[iter][i] - theta[iter-1][i]) 
	           /(fabs(theta[iter-1][i]) + REAL(epsilon1)[0]) ;
	R_rsort(diff, nB+2) ;
        if (diff[nB+1]<=REAL(epsilon2)[0]){
	  if (INTEGER(trace)[0]==1) 
	    Rprintf("\nAlgorithm has converged. \n") ;
	  conv = 0 ;
	  break ;
	}
        
        /** determine if need to increase sample size **/
	if (INTEGER(fixed_size)[0]==0 &&
            iter <INTEGER(sample_iter)[0] &&
            da->nS < INTEGER(max_size)[0]){
          hessGlmEst(theta[iter], da, VL, HQ, VQ) ;
          // approx inverse of var of theta
          // since beta varies little, it's insensible to inverse the whole near singular matrix
          // instead, I just consider the submatrix corresp to phi and p
	  for (i=0;i<2;i++){	    
	    for (j=0;j<2;j++){
	      HQ2[i][j] = HQ[i+nB][j+nB]  ;		    
	      VQ2[i][j] = VQ[i+nB][j+nB]  ;		    
	    }
	  }
          smat_inverse(VQ2, mc_var_inv, 2) ;
	  smat_multiply(HQ2, mc_var_inv, mc_var_inv2, 2) ;
	  smat_multiply(mc_var_inv2, HQ2, mc_var_inv, 2) ;	  
	  for (i=nB;i<nB+2;i++)
	    diff[i] = theta[iter-1][i] - theta[iter][i] ;
          T =0 ;
	  for (i=0;i<2;i++){
	    for (j=0;j<2;j++)
	      T += diff[nB+j] * mc_var_inv[j][i] * diff[nB+i]  ;		    
	  }	 
	  if (T<=C)
              nS = da->nS = imin2(da->nS+fprec(da->nS/REAL(ck)[0],0),
                                INTEGER(max_size)[0]) ;
	}
        
	if (INTEGER(trace)[0]==1 && INTEGER(fixed_size)[0]==0)
	    Rprintf(" sample size: %d", da->nS) ;    	  
	Rprintf("\n") ;
    }

    // fit glm after convergence
    tw_glm(da->X, da->Y, da->offset, da->weights, da->beta,
           da->p, da->link_power, eps_glm, nO, nB, iter_glm,
           xw, beta_old, zw, qraux, work, qty, &pivot);
    cpglm_fitted(da->beta, da, da->mu, da->eta);

    if (diff[nB+1]>REAL(epsilon2)[0]){
        warning("Algorithm did not converge: iteration limit reached.");
	conv = 1 ;
    }

    /********************************************
    ** return result
    1. coefficients(beta) 2. residuals 3. fitted.values
    4. linear.predictors  5. iter  6. weights(working) 
    7. converged  8. phi  9. p  10. theta  11. theta.all  
    12. hess  13. final.size  14. imp.weights 
     *********************************************/ 
    SEXP ans, ans1, ans2, ans3, ans4, ans5, ans6, ans7,
      ans8, ans9, ans10, ans11, ans12, ans13, ans14, names;
     
    PROTECT(ans=allocVector(VECSXP,14)) ;
    PROTECT(ans1=allocVector(REALSXP,nB));
    PROTECT(ans2=allocVector(REALSXP,nO));
    PROTECT(ans3=allocVector(REALSXP,nO));
    PROTECT(ans4=allocVector(REALSXP,nO));
    PROTECT(ans5=allocVector(INTSXP,1));
    PROTECT(ans6=allocVector(REALSXP,nO));
    PROTECT(ans7=allocVector(INTSXP,1));
    PROTECT(ans8=allocVector(REALSXP,1));
    PROTECT(ans9=allocVector(REALSXP,1));
    PROTECT(ans10=allocVector(REALSXP,nB+2));
    PROTECT(ans11=allocMatrix(REALSXP, iter+1, nB+2));
    PROTECT(ans12=allocMatrix(REALSXP, nB+2, nB+2));
    PROTECT(ans13=allocVector(INTSXP,1));
    PROTECT(ans14=allocVector(REALSXP,nS));

    for (i=0;i<nB;i++)
      REAL(ans1)[i] = da->beta[i] ;
    for (i=0;i<nO;i++){
      REAL(ans2)[i] = (da->Y[i] - da->mu[i])/
            mu_eta(da->eta[i], da->link_power) ;
      REAL(ans3)[i] = da->mu[i] ;
      REAL(ans4)[i] = da->eta[i] ;  
      REAL(ans6)[i] = sqrt(da->weights[i]/varFun(da->mu[i],da->p)) *
	        mu_eta(da->eta[i], da->link_power);
    }
    INTEGER(ans5)[0] = iter ;
    INTEGER(ans7)[0] = conv ;
    REAL(ans8)[0] = da-> phi ;
    REAL(ans9)[0] = da-> p ;

    for (j=0;j<nB+2;j++){
      REAL(ans10)[j] = theta[iter][j] ;
      for (i=0;i<iter+1;i++)
	REAL(ans11)[i+(iter+1)*j]= theta[i][j] ;
    }

    /********************************************
     ** compute information matrix
     *********************************************/ 
    if (INTEGER(trace)[0]==1) 
      Rprintf("Computing covariance matrix...\n") ;
    hessGlmEst(theta[iter], da, VL, HQ, VQ);        
    for (i=0;i<nB+2;i++){
      for (j=0;j<nB+2;j++)
	REAL(ans12)[i+(nB+2)*j] =  -VL[i][j] - HQ[i][j] ;
    }

    INTEGER(ans13)[0] = da->nS ;
    for (i=0;i<da->nS;i++)
      REAL(ans14)[i] = da->iw[i] ;

    SET_VECTOR_ELT(ans, 0, ans1);
    SET_VECTOR_ELT(ans, 1, ans2);
    SET_VECTOR_ELT(ans, 2, ans3);
    SET_VECTOR_ELT(ans, 3, ans4);
    SET_VECTOR_ELT(ans, 4, ans5);
    SET_VECTOR_ELT(ans, 5, ans6);
    SET_VECTOR_ELT(ans, 6, ans7);
    SET_VECTOR_ELT(ans, 7, ans8);
    SET_VECTOR_ELT(ans, 8, ans9);
    SET_VECTOR_ELT(ans, 9, ans10);
    SET_VECTOR_ELT(ans, 10, ans11);
    SET_VECTOR_ELT(ans, 11, ans12);
    SET_VECTOR_ELT(ans, 12, ans13);
    SET_VECTOR_ELT(ans, 13, ans14);

    // set names attributes
    PROTECT(names = allocVector(STRSXP, 14));
    SET_STRING_ELT(names, 0, mkChar("coefficients"));
    SET_STRING_ELT(names, 1, mkChar("residuals"));
    SET_STRING_ELT(names, 2, mkChar("fitted.values"));
    SET_STRING_ELT(names, 3, mkChar("linear.predictors"));
    SET_STRING_ELT(names, 4, mkChar("iter"));
    SET_STRING_ELT(names, 5, mkChar("weights"));
    SET_STRING_ELT(names, 6, mkChar("converged"));
    SET_STRING_ELT(names, 7, mkChar("phi"));
    SET_STRING_ELT(names, 8, mkChar("p"));
    SET_STRING_ELT(names, 9, mkChar("theta"));
    SET_STRING_ELT(names, 10, mkChar("theta.all"));
    SET_STRING_ELT(names, 11, mkChar("hess"));
    SET_STRING_ELT(names, 12, mkChar("final.size"));
    SET_STRING_ELT(names, 13, mkChar("imp.weights"));
    setAttrib(ans, R_NamesSymbol, names);
    UNPROTECT(16) ;
    
    return ans ;

}
