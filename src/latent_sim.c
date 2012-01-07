/************************************************/
/*   Rejection sampling for latent variable     */
/************************************************/

/**
 * @file latent_sim.c
 * @brief Function for implementing rejection sampling
 * of the latent variables as need in Bayesian updating.   
 * @author Wayne Zhang                            
*/


#include "cplm.h"
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
 * @param da a SEXP struct
 *
 * @return log difference
 */ 
double cplm_rlatT_diff (double par,  SEXP da){
    int kk = K_ELT(da)[0], *ygt0 = YPO_ELT(da) ;
    int k = ygt0[kk] ;
    double *Y = Y_ELT(da),  *wts =PWT_ELT(da),
        p = P_ELT(da)[0], phi = PHI_ELT(da)[0], lambda = LAM_ELT(da)[0];
    double phi2 = phi / wts[k] ;
    double ldiff = cplm_post_latT(par,Y[k],phi2, p)
        - cplm_dtpois(par,lambda) ;
    return ldiff ;
}

/**
 * Rejections sampling of the i_th latent variable T_i  
 * used in the cplm library
 *
 * @param da a list object
 *
 */
void cplm_rlatT_reject (SEXP da){
    int kk = K_ELT(da)[0], *ygt0 = YPO_ELT(da), *simT = SIMT_ELT(da) ;
    int k = ygt0[kk], xtemp ;
    double *Y = Y_ELT(da),  *wts =PWT_ELT(da),
        p = P_ELT(da)[0], phi = PHI_ELT(da)[0], *lambda = LAM_ELT(da);
    double par=2.0, val, pb,  phi2, mx;
    double c, b, p1= p -1, p2= 2 - p ;
    // adjust for prior weights
    phi2 = phi / wts[k] ;        
    *lambda = cplm_lambda_tpois(Y[k], phi2, p) ;

    // find t that results in the largest difference
    // between target and proposal
    c = p2/p1 * (log(Y[k]) - log(p1)) - log(phi2)/p1
        -log(p2)- log(*lambda) ;
    b = p2 / p1 ;
    if (c <0 ) 
        val = cplm_rlatT_diff(1.0, da) ;
    else {
        mx = exp((c-b*log(b))/b) ; //approximate mode
        if (mx <2)
            val = fmax2(fmax2(cplm_rlatT_diff(1.0, da),
                              cplm_rlatT_diff(2.0, da)),
                        cplm_rlatT_diff(3.0, da));
        else{
            par = ceil(mx) ;
            val = fmax2(cplm_rlatT_diff(par, da),
                        cplm_rlatT_diff(par-1.0, da));
        }
    }
    while (1) {
        xtemp =cplm_rtpois(*lambda);
        pb = exp(cplm_rlatT_diff(xtemp,da)-val);
        if (runif(0,1)<=pb){
            simT[kk] = xtemp ;
            break ;
        }
    }
}



/**
 * compute joint likelihood for one Single realization 
 * of simT (column vector)  given mu, phi and p
 *
 * @param da a SEXP struct
 *
 * @return log likelihood for one vector of T
 *
 */
double cplm_llik_lat(SEXP da){
    int *dm = DIMS_ELT(da), *ygt0 = YPO_ELT(da),*simT = SIMT_ELT(da) ;
    int nO = dm[nO_POS], nP = dm[nP_POS] ;
    double *Y = Y_ELT(da), *mu = MU_ELT(da), *wts =PWT_ELT(da),
    p = P_ELT(da)[0], phi =PHI_ELT(da)[0];    
    int i, k ;
    double lp=0, p2=2-p, p1=p-1;
    for (i=0; i<nO; i++)
        lp += pow(mu[i],p2) * wts[i];
    lp /= (- phi*p2) ;
    for (i=0; i<nP; i++){
	k = ygt0[i] ;
        lp += - Y[k]*pow(mu[k],-p1)*wts[k] /(phi*p1) +
            cplm_post_latT(simT[i],Y[k],phi/wts[k],p);
    }
    return lp ;
}        
