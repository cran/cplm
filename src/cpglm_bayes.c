/************************************************************/
/*   Function for the Markov Chain Monte Carlo algorithm    */
/*    in the Compound Poisson Generalized Linear Model      */
/*              Author:  Wayne Zhang                        */
/*            actuary_zhang@hotmail.com                     */
/************************************************************/

#include "cplm.h"

// struct used in cpglm_bayes
typedef struct {
  da_parm *dap ;          // struct to store data and parameters
  int *simTvec ;          // vector of simulated latent variables
  double *mu ;            // mean vector
  double *eta ;           // linear predictors
  double *pbeta_mean ;    // vector of prior means for beta
  double *pbeta_var ;	  // vector of prior variance for beta	
}  bcpglm_str;

/************************************************/
/*   Function to compute full conditionals      */  
/************************************************/

// posterior log density of the index parameter p
double bcpglm_post_p(double x, void *data){
  bcpglm_str *da = data ;
  double ld = cplm_llikS(da->mu, da->dap->phi, x,
			 da->simTvec, da->dap) ;
  return ld ;
}

/*
// posterior log density of betak - naive Gibbs update
double bcpglm_post_beta(double x, void *data){
  bcpglm_str *da = data ;
  da_parm *dap = da->dap ;
  int k = dap->k, i, kk ;
  double beta_old = dap->beta[k], 
    ld=0, p2=2-dap->p, p1=dap->p-1 ;
  dap->beta[k] = x ;
  // update mu
  cpglm_fitted(dap->beta, dap, da->mu, da->eta) ;
  // loglikelihood from data
  for (i=0; i<dap->nO; i++)
    ld += pow(da->mu[i],p2) * dap->weights[i];
  ld /= (- dap->phi*p2) ;
  for (i=0; i<dap->nP; i++){
    kk = dap->ygt0[i] ;
    ld += - dap->Y[kk]*pow(da->mu[kk],-p1)*
      dap->weights[kk] /(dap->phi*p1);
  }
  // prior info
  ld += -0.5*(x-da->pbeta_mean[k])*(x-da->pbeta_mean[k])/da->pbeta_var[k] ;
  // restore betak
  dap->beta[k] = beta_old ;
  return ld ;
}
*/

/* Not used as we use direct sampling from inverse Gamma for phi
// posterior log density of the dispersion parameter phi
double bcpglm_post_phi(double x, void *data){
  bcpglm_str *da = data ;
  da_parm *dap = da->dap ;
  int i, kk ;
  double ld=0, p2=2-dap->p, p1=dap->p-1 ;
  for (i=0; i<dap->nO; i++)
    ld += pow(da->mu[i],p2) * dap->weights[i];
  ld /= (- x*p2) ;
  for (i=0; i<dap->nP; i++){
    kk = dap->ygt0[i] ;
    ld += - dap->Y[kk]*pow(da->mu[kk],-p1)*dap->weights[kk] /(x*p1) 
      - log(x)*da->simTvec[i]/p1;
  }
  return ld ;
}
*/

// posterior log density of of the vector of beta
double bcpglm_post_beta_vec(double *x, void *data){
  bcpglm_str *da = data ;
  da_parm *dap = da->dap ;
  int i, kk ;
  double ld=0, p2=2-dap->p, p1=dap->p-1 ;
  // update mu
  cpglm_fitted(x, dap, da->mu, da->eta) ;
  // loglikelihood from data
  for (i=0; i<dap->nO; i++)
    ld += pow(da->mu[i],p2) * dap->weights[i];
  ld /= (- dap->phi*p2) ;
  for (i=0; i<dap->nP; i++){
    kk = dap->ygt0[i] ;
    ld += - dap->Y[kk]*pow(da->mu[kk],-p1)*
      dap->weights[kk] /(dap->phi*p1);
  }
  // prior info
  for (i=0;i<dap->nB;i++)
    ld += -0.5*(x[i]-da->pbeta_mean[i])*
      (x[i]-da->pbeta_mean[i])/da->pbeta_var[i] ;
  return ld ;
}

/************************************************/
/*   Multivariate Metropolis update with        */  
/*   random walk multivariate Normal proposal   */  
/************************************************/
/*
  simulation of multivariate normal
  d: dimension
  m: mean vector
  v: positive-definite covarince matrix
  s: vector to store simulated values
*/
// FIXME: does it handle d==1?
void rmvnorm(int d, double *m, double **v, double *s){
  int i, j, info, incx=1;
  double *lv = Calloc(d*d, double) ;
  char c1='L', c2='N' ;
  GetRNGstate() ;
  // simulate d univariate normal r.v.s
  for (i=0;i<d;i++)
    s[i] = rnorm(0,1) ;
  PutRNGstate() ;    
  // cholesky factor of v
  for (i=0;i<d;i++){
    for (j=0;j<d;j++)
      lv[i+j*d] = v[i][j] ;
  }    
  F77_CALL(dpotrf)(&c1,&d,lv,&d,&info) ;
  if (info!=0) 
    error("Error %d in Cholesky decomposition.", info) ;        
  // scale and shift univariate normal r.v.s
  F77_CALL(dtrmv)(&c1,&c2,&c2,&d,lv,&d,s,&incx) ;
  for (i=0;i<d;i++)
    s[i] += m[i] ;    
  Free(lv) ;    
}

/* random walk metropolis sampling for a vector of parameters 
   of length d using multivariate normal proposal
   sn: simulated new vector 
   myfunc: user specified function to compute log posterior 
   data: struct used in myfunc
   return value is a 0-1 integer: 0 means not accept
*/
int metrop_mvnorm_rw(int d, double *m, double **v, double *sn, 
		     double (*myfunc)(double *x, void *data), 
		     void *data){
  double A ;
  rmvnorm(d,m, v, sn) ;
  // determine if accept the sample
  A = exp(myfunc(sn,data)-myfunc(m,data) ) ;
  if (A<1 && runif(0,1)>=A)
    return 0 ;
  else 
    return 1 ;
}

/************************************************/
/*     Main function to fit compound Poisson    */
/*     GLM using Monte Carlo Markov Chains      */
/************************************************/

SEXP bcpglm_gibbs_lat (SEXP X,             // design matrix
		       SEXP Y,             // response vector
		       SEXP ygt0,          // row index of positive values 
		       SEXP offset,        // offset
		       SEXP weights,       // prior weights
                       SEXP beta,          // initial model coefficients 
		       SEXP phi,           // initial dispersion parameter
		       SEXP p,             // initial index parameter
		       SEXP link_power,    // power of link function, as in tweedie
                       SEXP pbeta_mean,    // vector of prior mean for beta
		       SEXP pbeta_var,     // vector of prior variance for beta
		       SEXP pphi_shape,    // prior shape parameter in Inverse-Gamma for phi
                       SEXP pphi_scale,    // prior scale parameter in Inverse-Gamma for phi
		       SEXP bound_p,       // vector of lower and upper bound of p
		       SEXP ebeta_var,     // cov matrix used in Metropolis update
                       SEXP niter,         // # iterations
		       SEXP nburnin,       // # burn in
		       SEXP nthin,         // thinning rate
		       SEXP nkeep,         // # simulations in final result
                       SEXP nreport)       // # report frequency
{

  // dimensions 
  int nO=LENGTH(Y),
    nP = LENGTH(ygt0),
    nB = LENGTH(beta),
    nit = INTEGER(niter)[0],
    nbn = INTEGER(nburnin)[0], 
    nth = INTEGER(nthin)[0],
    nS = INTEGER(nkeep)[0],
    nR = INTEGER(nreport)[0];    
  int iter, i, j, ns, kk ;
  int err, ninit = 4, acc, accept=0;  
  double xl_p = REAL(bound_p)[0], xr_p =REAL(bound_p)[1];
  double xtemp,  **sims;
  double shape = REAL(pphi_shape)[0], scale = REAL(pphi_scale)[0] ;
  double *beta_sim, **beta_norm_var;
  double  p1, p2, alp, gamm;
  // allocate memory for data struct
  bcpglm_str *da = (bcpglm_str *) R_alloc(1,sizeof(bcpglm_str)) ;
  da_parm *dap = (da_parm *) R_alloc(1,sizeof(da_parm)) ;
  da->dap = dap ;
  dap->X = dmatrix(nO,nB) ;
  dap->Y = dvect(nO) ;
  dap->ygt0 = ivect(nP) ;
  dap->offset = dvect(nO) ;
  dap->weights = dvect(nO) ;
  dap->beta= dvect(nB) ;
  da->mu = dvect(nO) ;
  da->eta = dvect(nO) ;
  da->simTvec = ivect(nP) ;
  da->pbeta_mean = dvect(nB) ;
  da->pbeta_var = dvect(nB) ;
  sims = dmatrix(nS,nB+2) ;
  beta_sim = dvect(nB) ;
  beta_norm_var = dmatrix(nB,nB) ;
     
  // fill in struct dap and da
  dap->nO = nO ;
  dap->nP = nP;    
  dap->nB = nB ;
  for (i=0;i<nP;i++)
    // ygt0 from R does not start from 0
    dap->ygt0[i]=INTEGER(ygt0)[i]-1 ;
  for (i=0;i<nO;i++){
    dap->Y[i]=REAL(Y)[i] ;
    dap->offset[i]=REAL(offset)[i] ;
    dap->weights[i]=REAL(weights)[i] ;
    for (j=0;j<nB;j++)
      dap->X[i][j] = REAL(X)[i+j*nO];
  }
  for (i=0;i<nB;i++){
    dap->beta[i] = REAL(beta)[i] ;
    da->pbeta_mean[i] = REAL(pbeta_mean)[i] ;
    da->pbeta_var[i] = REAL(pbeta_var)[i] ;
  }
  dap->phi = REAL(phi)[0] ;
  dap->p = REAL(p)[0];
  dap->link_power = REAL(link_power)[0] ;
  for (i=0;i<nB;i++){
    for (j=0;j<nB;j++)
      beta_norm_var[i][j] = REAL(ebeta_var)[i+j*nB];
  }    
  cpglm_fitted(dap->beta, dap, da->mu, da->eta);

  // MCMC iterations
  GetRNGstate() ;
  for (iter=0;iter<nit;iter++){
    if (nR>0 && (iter+1)%nR==0)
      Rprintf("Iteration: %d \n ", iter+1) ;
    R_CheckUserInterrupt() ;
        
    // update latent variable T using rejection sampling
    for (i=0;i<nP;i++){    
      da->dap->k=i ;
      cplm_rlatT_reject(1, &(da->simTvec[i]), da->dap) ;
    }        	
    R_CheckUserInterrupt() ;

    // update p using ARS (? need metropolis?)
    err = arms_simple(ninit,&xl_p,&xr_p,bcpglm_post_p,
		      (void *) da,1,&(dap->p),&xtemp);
    if (err>0) 
      error("Sample error %d in updating p\n",err) ;
    else 
      dap->p = xtemp ;    
    R_CheckUserInterrupt() ;
        
    // update mean parameters
    /*
    // naive gibbs update using ARS
    for (i=0;i<nB;i++){
    da->dap->k=i;
    xl_beta = dap->beta[i] - beta_fctr*fabs(dap->beta[i]) ;
    xr_beta = dap->beta[i] + beta_fctr*fabs(dap->beta[i]) ;
    err = arms_simple(ninit,&xl_beta,&xr_beta,bcpglm_post_beta,
    (void *) da, 0,&(dap->beta[i]),&xtemp);  
    if (err>0) 
    error("Sample error %d in updating beta\n",err) ;
    else 
    dap->beta[i] = xtemp ;
    }
    */
    //Metropolis-Hasting block update

    acc = metrop_mvnorm_rw(nB, dap->beta, beta_norm_var,
			   beta_sim, bcpglm_post_beta_vec, da) ;
    if (acc>0) {
      for (i=0;i<nB;i++)
	dap->beta[i] = beta_sim[i] ;
    }
    accept += acc ;
    if (nR>0 && (iter+1)%nR==0)
      Rprintf("Acceptance rate in Metropolis update: %4.2f%%\n",
	      accept*1.0/(iter+1)*100);                
    cpglm_fitted(dap->beta, dap, da->mu, da->eta) ;    
    R_CheckUserInterrupt() ;

    // direct simulation of phi from inverse-Gamma    
    p1 = dap->p -1 ;
    p2 = 2- dap->p ;
    alp = icumsum(da->simTvec, nP)/p1 + shape ;
    gamm =0 ;
    for (i=0; i<nO; i++)
      gamm += pow(da->mu[i],p2) * dap->weights[i];
    gamm /= p2 ;
    for (i=0; i<nP; i++){
      kk = dap->ygt0[i] ;
      gamm += dap->Y[kk]*pow(da->mu[kk],-p1)*dap->weights[kk] /p1 ;
    }
    gamm = 1/(gamm + scale) ; 
    dap->phi = 1/rgamma(alp,gamm);
    
    // store results 
    if (iter>=nbn &&  (iter+1-nbn)%nth==0 ){
      ns = (iter+1-nbn)/nth -1;   
      for (i=0;i<nB;i++)
	sims[ns][i] = dap->beta[i] ;
      sims[ns][nB] = dap->phi ;
      sims[ns][nB+1] = dap->p ;
      
    }
    
  }
  PutRNGstate() ;
    
  //return result    
  SEXP ans ;     
  PROTECT(ans=allocMatrix(REALSXP, nS, nB+2));
  for (j=0;j<nB+2;j++){
    for (i=0;i<nS;i++)
      REAL(ans)[i+nS*j]= sims[i][j] ;
  }
  UNPROTECT(1) ;
  return ans ;
    
}

