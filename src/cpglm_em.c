/************************************************************/
/*   Function to implement the Monte Carlo EM algorithm     */
/*    in the Compound Poisson Generalized Linear Model      */
/*              Author:  Wayne Zhang                        */
/*            actuary_zhang@hotmail.com                     */
/************************************************************/

/**
 * @file cpglm_em.c
 * @brief Function for implementing the MCEM algorithm
 * in the Compound Poisson Generalized Linear Model            
 * @author Wayne Zhang                         
 */

#include "cplm.h"

// struct to store info related to latent variables 
typedef struct {
    int nS ;               /**< number of simulations */
    int **simT ;           /**< 2d array of simulated latent Poisson variable */
    int *simTvec ;         /**< vector of one realization of the Poisson variable */
    double *et ;           /**< expected value of latent Poisson E(simT) */
    double *iw ;           /**< importance weights  */
    double sw ;            /**< sum of importance weights */
}  lat_bag;

// struct used in cpglm
typedef struct {
    da_parm *dap ;        /**< struct to store data and parameters */
    lat_bag *lat ;        /**< struct for latent variables */
    double *mu ;          /**< mean vector */
    double *eta ;         /**< linear predictor */
    double *mu1p ;        /**< mu^(1-p) */ 
    double *mueta ;	  /**< d(mu)/d(eta)*/	
}  cpglm_str;

/** relative tolerance in tweedie GLM */
#define EPS_GLM 1.0E-8


/************************************************/
/*   Rejection sampling for latent variable     */
/************************************************/

/** 
 * Log posterior density of the latent variable T_i
 *
 * @param x a point at which the density to be computed
 * @param y the positive response value corresponding to the t_i
 * @param phi scale parameter (if weights supplied, phi/weights)
 * @param p index parameter
 *
 * @return log density of T_i
 */

double cplm_post_latT(double x, double y, double phi, double p){
    double p1=p-1, p2=2-p, ans;
    ans = x*(p2/p1 * (log(y) - log(p1)) - log(phi)/p1 -log(p2)) -
        lgamma(x*p2/p1) - lgamma(x+1) ;
    return ans ;
}

/**
 * function to compute the mean of the proposal 0-truncated poisson 
 * by finding the approximate mode
 *
 * @param y obseved value from the tweedie distribution
 * @param phi scale parameter (if weights supplied, phi/weights)
 * @param p index parameter
 *
 * @return approximate mode of tweedie
 */

double cplm_lambda_tpois(double y, double phi, double p){
    double lamdba ;
    lamdba = ceil(pow(y,2-p)/((2-p)*phi));
    return lamdba ;
}

/**
 * function to simulate from 0 truncated poisson
 *
 * @param lambda mean of the untrucated poisson
 *
 * @return simulated value from 0-truncated poisson
 */
int cplm_rtpois(double lambda ){
    int sim = qpois(runif(dpois(0, lambda,0), 1), lambda, 1,0);
    return sim ;
}

/**
 * function to compute the density of 0-truncated poisson on log scale
 *
 * @param x the point at which the density is to be evaluated
 * @param lambda  mean of the original poisson
 *
 * @return log density
 */    
double cplm_dtpois(double x, double lambda){
    double ld ;
    ld = x*log(lambda)- lambda -lgamma(x+1)-log(1-exp(-lambda)) ;
    return ld ;
}

/**
 * compute difference between target and proposal (unnormalized)
 * used to locate the point where max diff occurs
 *
 * @param par the point at which to evaluate the difference
 * @param dap a da_parm struct
 *
 * @return log difference
 */ 
double cplm_rlatT_diff (double par,  da_parm *dap){
    int k = dap->ygt0[dap->k] ;
    double phi2 = dap->phi / dap->weights[k] ;
    double ldiff = cplm_post_latT(par,dap->Y[k],phi2, dap->p)
        - cplm_dtpois(par,dap->lambda) ;
    return ldiff ;
}

/**
 * Rejections sampling of the i_th latent variable T_i  
 * used in the cplm library
 *
 * @param nS number of simulated values
 * @param ans a vector of nS simulated values
 * @param dap a da_parm struct
 *
 */
void cplm_rlatT_reject (int nS, int *ans, da_parm *dap){
    int i, xtemp,  k=dap->ygt0[dap->k];
    double par=2.0, val, pb,  phi2, mx;
    double c, b, p1= dap->p -1, p2= 2 - dap->p ;
    // adjust for prior weights
    phi2 = dap->phi / dap->weights[k] ;        
    dap->lambda = cplm_lambda_tpois(dap->Y[k], phi2, dap->p) ;

    // find t that results in the largest difference
    // between target and proposal
    c = p2/p1 * (log(dap->Y[k]) - log(p1)) - log(phi2)/p1
        -log(p2)- log(dap->lambda) ;
    b = p2 / p1 ;
    if (c <0 ) 
        val = cplm_rlatT_diff(1.0, dap) ;
    else {
        mx = exp((c-b*log(b))/b) ; //approximate mode
        if (mx <2)
            val = fmax2(fmax2(cplm_rlatT_diff(1.0, dap),
                              cplm_rlatT_diff(2.0, dap)),
                        cplm_rlatT_diff(3.0, dap));
        else{
            par = ceil(mx) ;
            val = fmax2(cplm_rlatT_diff(par, dap),
                        cplm_rlatT_diff(par-1.0, dap));
        }
    }
    for (i=0;i<nS;i++){
        while (1) {
            xtemp =cplm_rtpois(dap->lambda);
            pb = exp(cplm_rlatT_diff(xtemp,dap)-val);
            if (runif(0,1)<=pb){
                ans[i] = xtemp ;
                break ;
            }
        }
    }
}

/**
 * Rejection sampling specific to cpglm_em
 *
 * @param ans a 2d array to store the simulatd results
 * @param da a cpglm_str struct
 *
 */
static void cpglm_rlatT_reject (int **ans, cpglm_str *da){
    int i ;
    GetRNGstate() ;
    for (i=0; i< da->dap->nP;i++){    
        da->dap->k=i ;
        cplm_rlatT_reject(da->lat->nS, ans[i], da->dap) ;
    }
    PutRNGstate() ;
}

/************************************************/
/*          GLM in the tweedie family           */
/************************************************/

/**
 * function to fit tweedie glm
 *
 * @param x design matrix in long format
 * @param y response variable
 * @param off  offsets
 * @param wts prior weights
 * @param beta starting values of beta
 * @param vp variance power
 * @param lp link power
 * @param n number of obs
 * @param p number of parameters
 *
 */
void cplm_tw_glm(double *x, double *y, double *off, double *wts, 
                 double *beta, double vp, double lp,int n, int p){
    
    int i,j, k, job, info, nf, nx=100, pivot ;
    double  dum,  eta, mu, mu_eta_val,z, w;
    double *beta_old=Calloc(n, double),
        *xw = Calloc(n*p, double),
        *zw = Calloc(n, double),
        *qraux = Calloc(p, double),
        *work = Calloc(p, double),
        *qty = Calloc(n, double) ;
        
    // iterative re-weighted least squares
    for (k=0;k<nx;k++){
        // compute eta
        for (i=0;i<n;i++){    
            eta = 0 ;
            for (j=0;j<p;j++)
                eta += x[i+j*n] * beta[j] ; 
            eta += off[i] ;
            mu = linkInv(eta,lp) ;
            mu_eta_val = mu_eta(eta,lp);
            // linearization to get new response z
            z = eta + (y[i] - mu) / mu_eta_val -off[i] ;
            // working weight in least square
            w = sqrt(wts[i]/varFun(mu,vp) ) * mu_eta_val;
            // xw and zw are the design matrix and response in weighted ls
            for (j=0;j<p;j++)
                xw[i+j*n] = x[i+j*n] * w ;
            zw[i] = z * w ;
        }
        // old beta
        Memcpy(beta_old, beta, p) ;
        // qr decomposition (stored in xw now)
        job = 0 ; // no pivoting
        F77_CALL(dqrdc)(xw, &n, &n, &p, qraux, &pivot, work, &job) ;
        // weighted least square fit
        job = 100; // only want to compute beta
                   // but qty is also computed automatically by dqrsl
        F77_CALL(dqrsl)(xw, &n, &n, &p, qraux, zw,
                        &dum, qty, beta, &dum, &dum, &job, &info) ;

        // check if non-finite value is generated
        nf=0 ; 
        for (i=0;i<p;i++)
            nf += R_FINITE(beta[i]) ;
        if (nf<p)
            error(_("non-finite coefficients are found \n"));
        // branch out if converged
        if (dist(beta, beta_old, p)/ norm(beta, p) < EPS_GLM)
            break ;    
    }
    
    if (dist(beta, beta_old, p)/ norm(beta,p) >= EPS_GLM)
        warning("GLM did not converge - max iteration numbers reached!\n");

    Free(beta_old);
    Free(xw);
    Free(zw);
    Free(qraux);
    Free(work);
    Free(qty);
}
    
/************************************************/
/*   Single, Expected joint loglikelihood       */
/*   Single, Expected gradient and hessian      */
/*     of the joint loglikelihood               */
/************************************************/

/**
 * compute joint likelihood for one Single realization 
 * of simT (column vector)  given mu, phi and p
 *
 * @param mu mean vector
 * @param phi scale parameter
 * @param p index parameter
 * @param simT one realization of the simulated latent variables 
 * @param dap a da_parm struct
 *
 * @return log likelihood for one vector of T
 *
 */
double cplm_llikS(double *mu, double phi, double p,
                  int *simT, da_parm *dap){
    int i, k ;
    double lp=0, p2=2-p, p1=p-1;
    for (i=0; i<dap->nO; i++)
        lp += pow(mu[i],p2) * dap->weights[i];
    lp /= (- phi*p2) ;
    for (i=0; i<dap->nP; i++){
	k = dap->ygt0[i] ;
        lp += - dap->Y[k]*pow(mu[k],-p1)*dap->weights[k] /(phi*p1) +
            cplm_post_latT(simT[i],dap->Y[k],phi/dap->weights[k],p);
    }
    return lp ;
}

/**
 * Single realizion likelihood in cpglm
 *
 * @param mu mean vector
 * @param phi scale parameter
 * @param p index parameter
 * @param da a cpglm_str struct
 *
 * @return log likelihood for one vector of T
 */
static double cpglm_llikS(double *mu, double phi,
                          double p, cpglm_str *da){
    da_parm *dap = da->dap ;    
    int i ;
    double lp ;
    for (i=0;i<dap->nP;i++)
        da->lat->simTvec[i] = da->lat->simT[i][da->dap->k] ;
    lp = cplm_llikS(mu, phi, p, da->lat->simTvec, dap) ;
    return lp ;    
}

/**
 * Expected joint likelihood in cpglm
 *
 * @param mu mean vector
 * @param phi scale parameter
 * @param p index parameter
 * @param da a cpglm_str struct
 *
 * @return expected log likelihood, averaged over all Ts
 */
static double cpglm_llikE(double *mu, double phi, double p, cpglm_str *da){
    da_parm *dap = da->dap ;
    lat_bag *lat = da->lat ;
    int i, j, k;
    double lp=0, elgt, p1=p-1, p2=2-p;

    // taking expectation of dJLogLikS
    for (i=0; i<dap->nO; i++)
        lp += pow(mu[i],p2) * dap->weights[i];
    lp /= - phi*p2 ;
    for (i=0; i<dap->nP; i++){
        k = dap->ygt0[i] ;
        elgt =0 ; 
        for (j=0; j<lat->nS;j++)
            elgt += lgammafn(lat->simT[i][j]*p2/p1)* lat->iw[j] ;         
        elgt /= lat->sw ;
        lp += - dap->Y[k]*pow(mu[k],-p1)*dap->weights[k]/(phi*p1)- elgt
            + lat->et[i]*(p2/p1*(log(dap->Y[k])-log(p1))-log(phi/dap->weights[k])/p1 - log(p2));
    }    
    return lp ;
}

/**
 * compute derivative for Single L given one simulation of T:
 * beta and phi is computed analytically and p numerically
 *
 * @param par vector that stores beta, phi and p in order
 * @param gr vector to store gradient
 * @param da a cpglm_str struct
 *
 */
static void cpglm_llik_gradS (double *par, double *gr, cpglm_str *da){
    da_parm *dap = da->dap ;
    lat_bag *lat = da->lat ;
    int i, j, k, nB=dap->nB, nO=dap->nO ;
    double phi=par[nB], p=par[nB+1], p2=2-p, p1=p-1;

    // derivative for beta
    for (j=0; j<nB;j++){
        gr[j] = 0 ;
        for (i=0;i<dap->nO;i++)
            gr[j] += - dap->weights[i]* da->mu1p[i]
                * da->mueta[i]*dap->X[i+j*nO];
        for (i=0;i<dap->nP;i++){
            k = dap->ygt0[i] ;
            gr[j] += dap->Y[k]* da->mu1p[k]* dap->weights[k]
                /da->mu[k]* da->mueta[k]*dap->X[k+j*nO] ;
        }
        gr[j] /= phi ;
    }

    // derivative for phi
    gr[nB]=0 ;
    for (i=0;i<dap->nO;i++)
        gr[nB] += dap->weights[i]*da->mu1p[i]*da->mu[i];
    gr[nB] /= (phi*phi*p2);
    for (i=0;i<dap->nP;i++) {
        k = dap->ygt0[i] ;
        gr[nB] += dap->weights[k]*dap->Y[k]*da->mu1p[k]/(p1*phi*phi)
            - lat->simT[i][dap->k]/(p1*phi) ;
    }

    // derivative for p, using numerical solution
    //FIXME: did not consider bound here
    double dp1, dp2 ;
    p += EPS ;
    dp1 = cpglm_llikS(da->mu,phi,p,da) ;
    p -= 2*EPS ;
    dp2 = cpglm_llikS(da->mu,phi,p,da) ;
    gr[nB+1]=(dp1-dp2)/(2*EPS);
}

/**
 * compute  derivative for E(L), MC average across all T's
 *
 * @param par vector that stores beta, phi and p in order
 * @param gr vector to store gradient
 * @param da a cpglm_str struct
 *
 */ 
static void cpglm_llik_gradE (double *par, double *gr,  cpglm_str *da){
    da_parm *dap = da->dap ;
    lat_bag *lat = da->lat ;
    int i, j, k, nB=dap->nB, nO=dap->nO ;
    double phi=par[nB], p=par[nB+1];
    double p2=2-p, p1=p-1;
    
    // derivative for beta
    for (j=0; j<nB;j++){
        gr[j] = 0 ;
        for (i=0;i<dap->nO;i++)
            gr[j] += - dap->weights[i] * da->mu1p[i]
                * da->mueta[i] * dap->X[i+j*nO];
        for (i=0;i<dap->nP;i++){
            k = dap->ygt0[i] ;
            gr[j] += dap->weights[k] * dap->Y[k]* da->mu1p[k]/ da-> mu[k]
                * da->mueta[k] * dap->X[k+j*nO] ;
        }
        gr[j] /= phi ;
    }

    // derivative for phi
    gr[nB]=0 ;
    for (i=0;i<dap->nO;i++)
        gr[nB] += dap->weights[i]*da->mu1p[i]*da->mu[i];
    gr[nB] /= phi*phi*p2 ;
    for (i=0;i<dap->nP;i++) {
        k = dap->ygt0[i] ;
        gr[nB] += dap->weights[k]*dap->Y[k]* da->mu1p[k]/(p1*phi*phi)
            - lat->et[i]/(p1*phi) ;
    }

    // derivative for p, using numerical solution
    // FIXME: did not consider bound here
    double dp1, dp2 ;
    p += EPS ;
    dp1 = cpglm_llikE(da->mu, phi, p, da) ;
    p -= 2*EPS ;
    dp2 = cpglm_llikE(da->mu, phi, p, da) ;
    gr[nB+1]=(dp1-dp2)/(2*EPS);
}


/**
 * update mu, eta, mu1p, and mueta in cpglm_str
 * @param par vector that stores beta, phi and p in order
 * @param da a cpglm_str struct
 *
 */
static void cpglm_update_mustat (double *par, cpglm_str *da){
    da_parm *dap= da->dap ;
    int i ;
    double p = par[dap->nB+1], *beta_old=dap->beta ;
    dap->beta = par ;
    cpglm_fitted(da->eta, da->mu, da->mueta, dap);
    dap->beta = beta_old ;
    for (i=0;i<dap->nO;i++)
        da->mu1p[i]= pow(da->mu[i],1-p) ;
}

/**
 * compute hessian numerically
 *
 * @param par vector that stores beta, phi and p in order
 * @param hess matrix to store hessian
 * @param da a cpglm_str struct
 * 
 */
static void cpglm_llik_hessE (double *par, double **hess, cpglm_str *da){
    int i, j, nB=da->dap->nB ;    
    double *df1 = dvect(nB+2);
    double *df2 = dvect(nB+2);

    for (i = 0; i < nB+2; i++) {
	par[i] +=  EPS;
        cpglm_update_mustat(par, da) ;
	cpglm_llik_gradE(par, df1, da);
	par[i] -=  2 * EPS;
        cpglm_update_mustat(par, da) ;        
        cpglm_llik_gradE(par, df2, da);
	for (j = 0; j < nB+2; j++)
	    hess[i][j] = (df1[j] - df2[j])/ (2*EPS) ;
	par[i] +=EPS;
    }
    cpglm_update_mustat(par, da) ; // restore the stats       
}

/************************************************/
/*   Conditional likelihood and gradient        */
/*   for the index parameter                    */
/************************************************/
/**
 *
 * function to compute the (expected) log posterior of p
 *
 * @param x value at which to compute the log density
 * @param da a cpglm_str struct
 *
 * @return log posterior density of p
 */
static double cpglm_post_p(double x, cpglm_str *da){
    return cpglm_llikE(da->mu, da->dap->phi, x, da) ;
}

/**
 * Compute the log posterior of p and its gradient, as required in lfgsbU
 *
 * @param n not used but required in lfgsbU
 * @param par vector that stores beta, phi and p in order
 * @param data a void struct
 *
 */
static double cpglm_post_p_opt(int n, double *par, void *data){
    cpglm_str *da =data ;    
    double ans= cpglm_post_p(*par, da) * SCALE;
    return ans ;
}

/**
 * compute numerical derivative
 *
 * @param n not used but required in lfgsbU
 * @param par vector that stores beta, phi and p in order
 * @param gr gradient vector
 * @param data a void struct
 */
static void cpglm_post_p_grad_opt (int n, double *par, double *gr,  void *data){
    cpglm_str *da =data ;
    double par1 = *par +EPS, par2 =*par-EPS ;
    // FIXME: did not consider bound here
    *gr=(cpglm_post_p(par1,da)-cpglm_post_p(par2,da))/(2*EPS) * SCALE;
}

/************************************************/
/*   Function to compute variance matrix        */  
/*        using Monte Carlo methods             */
/************************************************/
/**
 * Function to compute variance matrix using Monte Carlo methods
 *
 * @param par all the parameters in model: beta, phi,p in order
 * @param da cpglm_str
 * @param VL variance of L
 * @param HQ hessian of E(L)
 * @param VQ approx variance of E(L)
 */
static void cpglm_hess (double *par,  cpglm_str *da,
                        double **VL, double **HQ, double **VQ){
    da_parm *dap = da->dap ;
    lat_bag *lat = da->lat ;
    int nS=lat->nS, nB = dap->nB ;
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
        dap->k=i ;
        cpglm_llik_gradS(par, grad, da);
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
    cpglm_llik_hessE(par, HQ, da);            
}



/************************************************/
/*   Utility functions used in cpglm_em         */  
/************************************************/
/**
 * Compute weights in importance sampling
 *
 * @param phi_new new value of phi
 * @param p_new new value of p
 * @param phi_old old value of phi
 * @param p_old old value of p
 * @param da a cpglm_str struct
 * @param iw vector to store importance weights
 * @param sw sum of importance weights
 */
static void cpglm_import_weight(double phi_new, double p_new, double phi_old,
                                double p_old, cpglm_str *da, double *iw, double *sw){
    da_parm *dap = da->dap ;
    lat_bag *lat = da->lat ;
    int nS = lat->nS, nP = dap->nP ;
    int i, j, k;
    double iw_max=-1.0E20 ;
    // compute weight
    for (j=0;j<nS;j++){
        iw[j] =0.0 ;
        for (i=0;i<nP;i++){
            k = dap->ygt0[i] ;
            iw[j] += cplm_post_latT(lat->simT[i][j],dap->Y[k], phi_new, p_new)
                - cplm_post_latT(lat->simT[i][j],dap->Y[k], phi_old, p_old);
        }
        iw_max = fmax2(iw_max, iw[j]) ;    
    }
    for (j=0;j<nS;j++)
        iw[j] = exp(iw[j]-iw_max) ; //subtract the max to avoid numerical issue
    *sw = dcumsum(iw, nS) ;
}

/**
 * initialize importance weight to unity
 *
 * @param nS number of simulations
 * @param iw importance weights
 * @param sw sum of iw
 */
static void unit_weight(int nS, double *iw, double *sw){
    int i ;
    for (i=0;i<nS;i++)
        iw[i] = 1.0 ;
    *sw = dcumsum(iw,nS) ;
}

/**
 * store parameter values to a matrix theta
 *
 * @param theta vector of optimized parameters
 * @param da a da_parms struct
 */
static void cpglm_store_theta(double *theta, da_parm *da){
    int i ;
    for (i=0;i< da->nB;i++)
        theta[i]= da->beta[i] ;
    theta[da->nB] = da->phi ;
    theta[da->nB+1] = da->p ;
}

/************************************************/
/*     Main function to fit compound Poisson    */
/*     GLM using Monte Carlo EM algorithm       */
/************************************************/

/**
 * Main function to fit compound Poisson GLM using Monte Carlo EM algorithm
 * @param X design matrix
 * @param Y response vector
 * @param ygt0  row index of positive values 
 * @param offset  offset
 * @param weights  prior weights
 * @param beta initial model coefficients 
 * @param phi initial dispersion parameter
 * @param p initial index parameter
 * @param link_power power of link function, as in tweedie
 * @param bound a vector giving the lower and upper bound of p
 * @param init_size initial sample size
 * @param sample_iter # iterations after which to use importance sampling
 * @param max_iter max # iterations and max sample size
 * @param epsilon1 stopping rule parameter
 * @param epsilon2 stopping rule parameter
 * @param ck stopping rule parameter
 * @param fixed_size whether sample size should be fixed
 * @param trace if fitting info should be printed
 * @param max_size max # samples 
 * @param beta_step integer, number of iterations to skip for updating beta
 *
 * @return a R list 
 */

SEXP cpglm_em (SEXP X,             // design matrix
	       SEXP Y,             // response vector
	       SEXP ygt0,          // row index of positive values 
	       SEXP offset,        // offset
	       SEXP weights,       // prior weights
	       SEXP beta,          // initial model coefficients 
	       SEXP phi,           // initial dispersion parameter
	       SEXP p,             // initial index parameter
	       SEXP link_power,    // power of link function, as in tweedie
	       SEXP bound,         // a vector giving the lower and upper bound of p
               SEXP init_size,     // initial sample size
	       SEXP sample_iter,   // # iterations after which to use importance sampling
	       SEXP max_iter,      // max # iterations and max sample size
	       SEXP epsilon1,      // stopping rule parameter
	       SEXP epsilon2,      // stopping rule parameter
	       SEXP ck,            // stopping rule parameter
	       SEXP fixed_size,    // whether sample size should be fixed
	       SEXP trace,         // if fitting info should be printed
	       SEXP max_size,      // max # samples 
	       SEXP beta_step)     // integer, number of iterations to skip for updating beta
{
    // dimensions 
    int nO=LENGTH(offset),
        nP = LENGTH(ygt0),
        nS=INTEGER(init_size)[0],
        nB=LENGTH(beta);
    int i, j, k, ii, conv, iter=0;    
    // em iteration variables 
    double **theta, *diff, *diff2, p1, p2, dpval, pnew;
    // var-cov related variables
    double  **VQ, **HQ, **VL;

    // allocate memory for data struct
    cpglm_str *da = (cpglm_str *) R_alloc(1,sizeof(cpglm_str)) ;
    da_parm *dap = (da_parm *) R_alloc(1,sizeof(da_parm)) ;
    lat_bag *lat = (lat_bag *) R_alloc(1,sizeof(lat_bag)) ;
    da->dap = dap ;
    da->lat = lat ;
    lat->et = dvect(nP) ;
    lat->simTvec = ivect(nP) ;
    da->eta = dvect(nO) ;
    da->mu = dvect(nO) ;
    da->mu1p = dvect(nO) ;
    da->mueta = dvect(nO) ;
 
    if (INTEGER(fixed_size)[0]==1){ 
        lat->simT = imatrix(nP, nS) ;
        lat->iw = dvect(nS) ;
    }
    else {// set column size to max possible to avoid resizing
        lat->simT = imatrix(nP, INTEGER(max_size)[0]) ; 
        lat->iw = dvect(INTEGER(max_size)[0]) ;
    }
    // initialize to unit weight
    unit_weight(nS, lat->iw, &(lat->sw));
 
    // for em iterations
    diff = dvect(nB+2) ;
    diff2 = dvect(INTEGER(max_iter)[0]) ;
    VL = dmatrix(nB+2, nB+2) ;
    HQ = dmatrix(nB+2, nB+2) ;
    VQ = dmatrix(nB+2, nB+2) ;
    
    // row dimension as max possible iterations
    theta = dmatrix(INTEGER(max_iter)[0],nB+2) ;
 
    // fill in struct dap, lat and da
    dap->nO = nO ;
    dap->nP = nP;    
    dap->nB = nB ;
    lat->nS = nS ;
    dap->Y = REAL(Y) ;
    dap->ygt0=INTEGER(ygt0);
    dap->offset=REAL(offset) ;
    dap->weights=REAL(weights) ;
    dap->X = REAL(X);
    dap->beta = REAL(beta) ;
    dap->phi = REAL(phi)[0] ;
    dap->p = REAL(p)[0];
    dap->link_power = REAL(link_power)[0] ;
    
    cpglm_fitted(da->eta, da->mu, da->mueta, dap);
 
    // store theta to begin with
    iter =0 ;
    cpglm_store_theta(theta[iter], dap) ;

    for (ii=0;ii<INTEGER(max_iter)[0];ii++){
        R_CheckUserInterrupt() ;
        iter++ ;
        /********************************************/
        /*************   E-step	   ******************/
        /********************************************/
        
	// Rejection sampling of T 
	if (iter <=INTEGER(sample_iter)[0]){	  
            cpglm_rlatT_reject(lat->simT,da) ;       
            if (INTEGER(fixed_size)[0]==0)
                unit_weight(lat->nS, lat->iw, &(lat->sw));                            
	}
	// Importance sampling of T 
	else {  
            cpglm_import_weight(dap->phi, dap->p, 
                                theta[INTEGER(sample_iter)[0]-1][nB],
                                theta[INTEGER(sample_iter)[0]-1][nB+1],
                                da, lat->iw, &(lat->sw)) ;
	}	  	  
	// compute weighted E(T_i)
	for (i=0;i<nP;i++)
            lat->et[i]=icumwsum(lat->simT[i],lat->iw,lat->nS)/lat->sw;

  	/********************************************/
        /****************   M-step    ***************/	
        /*********************************************/
                  
        // update p       
	pnew = dap->p ;
	lbfgsbU(&pnew, REAL(bound)[0], REAL(bound)[1], &dpval,
		cpglm_post_p_opt, cpglm_post_p_grad_opt,
                (void *)da,  &conv) ;
	// if converged, set to new value; o/w, retain
	if (conv == 0) 
            dap->p = pnew ;
        R_CheckUserInterrupt() ;
        
	// GLM update of beta 
        // do this every (beta_step) iterations as beta varies slowly        
        if ((iter-1)%INTEGER(beta_step)[0] ==0){
            cplm_tw_glm(dap->X, dap->Y, dap->offset, dap->weights, dap->beta,
                        dap->p, dap->link_power, nO, nB);
            // compute fitted values and related quantities 
            cpglm_fitted(da->eta, da->mu, da->mueta, dap);            
        }
        for (i=0;i<nO;i++)
            da->mu1p[i]= pow(da->mu[i],1-dap->p) ;
        
        // update phi         
	p1=dap->p-1, p2=2-dap->p;
        dap->phi = 0;
	for (i=0;i<nO;i++)
            dap->phi += pow(da->mu[i],p2)*dap->weights[i] ;
	dap->phi *= p1/p2 ;
	for (i=0;i<nP;i++){
            k = dap->ygt0[i] ;
            dap->phi += dap->Y[k]*dap->weights[k]*pow(da->mu[k],-p1);
	}
	dap->phi /= dcumsum(lat->et,nP) ;
        
        // store results 
        cpglm_store_theta(theta[iter], dap) ;
                
        // print out iteration info if necessary  
        if (INTEGER(trace)[0]==1){
            if (iter==INTEGER(sample_iter)[0]+1)
                Rprintf("Importance sampling begins...\n") ;
            Rprintf("Iteration: %d, phi: %6.3f, p: %6.4f",iter,dap->phi, dap->p) ;        
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

        // determine if need to increase sample size 
	if (INTEGER(fixed_size)[0]==0 &&
            iter <INTEGER(sample_iter)[0] &&
            lat->nS < INTEGER(max_size)[0]){
            
            for (i=nB;i<nB+2;i++)
                diff2[iter] += (theta[iter-1][i] - theta[iter][i])*
                    (theta[iter-1][i] - theta[iter][i]);
            diff2[iter] = sqrt(diff2[iter]) ;
            if (iter>=2){
                if (runif(0,diff2[iter-1]/diff2[iter])>0.5)
                    nS = lat->nS = imin2(lat->nS+fprec(lat->nS/REAL(ck)[0],0),
                                         INTEGER(max_size)[0]) ;
            }
	}
        
	if (INTEGER(trace)[0]==1 && INTEGER(fixed_size)[0]==0){
	    Rprintf(" sample size: %d", lat->nS) ;    	  
	    Rprintf("\n") ;
      }
        
    }
    
    // fit glm after convergence
    cplm_tw_glm(dap->X, dap->Y, dap->offset, dap->weights, dap->beta,
                dap->p, dap->link_power, nO, nB);
    cpglm_fitted(da->eta, da->mu, da->mueta,dap );
    
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
    PROTECT(ans14=allocVector(REALSXP,lat->nS));

    for (i=0;i<nB;i++)
        REAL(ans1)[i] = dap->beta[i] ;
    for (i=0;i<nO;i++){
        REAL(ans2)[i] = (dap->Y[i] - da->mu[i])/
            mu_eta(da->eta[i], dap->link_power) ;
        REAL(ans3)[i] = da->mu[i] ;
        REAL(ans4)[i] = da->eta[i] ;  
        REAL(ans6)[i] = sqrt(dap->weights[i]/varFun(da->mu[i],dap->p)) *
            mu_eta(da->eta[i], dap->link_power);
    }
    INTEGER(ans5)[0] = iter ;
    INTEGER(ans7)[0] = conv ;
    REAL(ans8)[0] = dap->phi ;
    REAL(ans9)[0] = dap->p ;

    for (j=0;j<nB+2;j++){
        REAL(ans10)[j] = theta[iter][j] ;
        for (i=0;i<iter+1;i++)
            REAL(ans11)[i+(iter+1)*j]= theta[i][j] ;
    }
    
    /*********************************************/
    /**      compute information matrix         **/
    /*********************************************/
    if (INTEGER(trace)[0]==1) 
        Rprintf("Computing covariance matrix...\n") ;
    cpglm_hess(theta[iter], da, VL, HQ, VQ);        
    for (i=0;i<nB+2;i++){
        for (j=0;j<nB+2;j++)
            REAL(ans12)[i+(nB+2)*j] =  -VL[i][j] - HQ[i][j] ;
    }

    INTEGER(ans13)[0] = lat->nS ;
    for (i=0;i<lat->nS;i++)
        REAL(ans14)[i] = lat->iw[i] ;

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
