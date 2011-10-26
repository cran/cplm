/************************************************************/
/*   Function for the Markov Chain Monte Carlo algorithm    */
/*    in the Compound Poisson Generalized Linear Model      */
/*              Author:  Wayne Zhang                        */
/*            actuary_zhang@hotmail.com                     */
/************************************************************/

/**
 * @file cpglm_bayes_lat.c
 * @brief Function for implementing the MCMC algorithm
 * in the Compound Poisson Generalized Linear Model            
 * @author Wayne Zhang                         
 */

#include "cplm.h"


/** struct used in cpglm_bayes */
typedef struct {
    da_parm *dap ;          /**< struct to store data and parameters  */
    int *simTvec ;          /**< vector of simulated latent variables */
    double *pbeta_mean ;    /**< vector of prior means for beta       */
    double *pbeta_var ;	    /**< vector of prior variance for beta    */
    double *mu ;            /**< mean vector */
    double *eta ;           /**< linear predictor */
    double bound_phi ;      /**< bound for phi */
    double *bound_p ;       /**< bound for p */    
    double mh_p_var ;       /**< proposal variance in M-H update for p */
    double mh_phi_var ;     /**< proposal variance in M-H update for phi */
    double *mh_beta_var ;   /**< proposal covariance matrix in M-H update for beta */
}  bcpglm_str;


/************************************************/
/*   Function to compute full conditionals      */  
/************************************************/

/**
 * posterior log density of the index parameter p
 *
 * @param x value of p at which the log density is to be calculated
 * @param data a void struct, cocerced to bcpglm_str internally
 *
 * @return log posterior density
 */
static double bcpglm_post_p_lat(double x, void *data){
    bcpglm_str *da = data ;
    double ld = cplm_llikS(da->mu, da->dap->phi, x,
                           da->simTvec, da->dap) ;
    return ld ;
}

/**
 * posterior log density of the dispersion parameter phi
 *
 * @param x value of phi at which the log density is to be calculated
 * @param data a void struct, cocerced to bcpglm_str internally
 *
 * @return log posterior density
 */
static double bcpglm_post_phi_lat(double x, void *data){
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
      - log(x)* da->simTvec[i]/p1;
  }
  return ld ;
}


/**
 * posterior log density of of the vector of beta
 *
 * @param x vector of values for beta
 * @param data void struct that is coerced to bcpglm_str
 *
 * @return log posterior density for beta
 */
double bcpglm_post_beta_lat(double *x,  void *data){
    bcpglm_str *da = data ;
    da_parm *dap = da->dap ;
    int i, kk, nO=dap->nO, nB=dap->nB ;
    double ld=0, p2=2-dap->p, p1=dap->p-1,
        *beta_old = dap->beta;
    // update mu
    dap->beta = x ;
    cpglm_fitted(da->eta, da->mu, (double *) NULL, dap);
    dap->beta = beta_old ;

    // loglikelihood from data
    for (i=0; i<nO; i++)
        ld += pow(da->mu[i],p2) * dap->weights[i];
    ld /= (- dap->phi*p2) ;
    for (i=0; i<dap->nP; i++){
        kk = dap->ygt0[i] ;
        ld += - dap->Y[kk]*pow(da->mu[kk],-p1)*
            dap->weights[kk] /(dap->phi*p1);
    }
    // prior info
    for (i=0;i<nB;i++)
        ld += -0.5*(x[i]-da->pbeta_mean[i])*
            (x[i]-da->pbeta_mean[i])/da->pbeta_var[i] ;
    return ld ;
}


/************************************************/
/*     Main function to fit compound Poisson    */
/*     GLM using Monte Carlo Markov Chains      */
/************************************************/

/**
 * MCMC simulation for compound Poisson GLM
 *
 * @param da a bcpglm_str struct
 * @param nR report interval
 * @param nit number of iterations
 * @param nbn number of burn-ins
 * @param nth thinning rate
 * @param sims a 2d array to store simulation results
 * @param acc_pct acceptance percentage 
 *
 */
static void bcpglm_mcmc_lat(bcpglm_str *da, int nR, int nit, int nbn,
                        int nth, double **sims, double *acc_pct){
    da_parm *dap = da->dap ;
    int nP = dap->nP, nB=dap->nB ; 
    double xtemp, *beta_sim ;
    double xl_p = da->bound_p[0], xr_p=da->bound_p[1],
        xr_phi = da->bound_phi, p_sd = sqrt(da->mh_p_var),
        phi_sd= sqrt(da->mh_phi_var) ;
    int  acc, accept[]={0,0,0}; 
    int i, j, iter, ns ;
    beta_sim = Calloc(nB, double) ;

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
        
        // M-H update of p using truncated normal
        acc = metrop_tnorm_rw(dap->p, p_sd, xl_p, xr_p, &xtemp, 
                                bcpglm_post_p_lat, (void *) da);
        dap->p = xtemp ;
        accept[0] += acc ;
        R_CheckUserInterrupt() ;
        
        //Metropolis-Hasting block update              
        acc = metrop_mvnorm_rw(nB, dap->beta, da->mh_beta_var,
                               beta_sim, bcpglm_post_beta_lat, (void *)da) ;
        Memcpy(dap->beta, beta_sim, nB) ;
        accept[1] += acc ;
        cpglm_fitted(da->eta, da->mu, (double*) NULL, dap) ;
        R_CheckUserInterrupt() ;
    

        // M-H update of phi using truncated normal
        acc = metrop_tnorm_rw(dap->phi, phi_sd, 0, xr_phi, &xtemp, 
                                bcpglm_post_phi_lat, (void *) da);
        dap->phi = xtemp ;
        accept[2] += acc ;
        R_CheckUserInterrupt() ;
        
        // print out acceptance rate if necessary
        if (nR>0 && (iter+1)%nR==0){
            Rprintf(_("Acceptance rate: beta(%4.2f%%), phi(%4.2f%%), p(%4.2f%%),\n"),
                    accept[1]*1.0/(iter+1)*100, accept[2]*1.0/(iter+1)*100,
                    accept[0]*1.0/(iter+1)*100 );
        }   
        // store results 
        if (iter>=nbn &&  (iter+1-nbn)%nth==0 ){
            ns = (iter+1-nbn)/nth -1;   
            for (j=0;j<nB;j++)
                sims[ns][j] = dap->beta[j];
            sims[ns][nB] = dap->phi  ;
            sims[ns][nB+1] = dap->p ;      
        } 
    }
    PutRNGstate() ;
    Free(beta_sim);
    // compute acceptance percentage
    for (i=0;i<3;i++)
        acc_pct[i] = accept[i]*1.0/nit ;
}


/**
 * implement MCMC for compound Poisson GLM using latent variables
 *
 * @param x a list object
 *
 * @return the simulated values 
 *
 */

SEXP bcpglm_gibbs_lat (SEXP x){
    // get dimensions
    int *dm = DIMS_ELT(x) ;
    int nO = dm[nO_POS],
        nP = dm[nP_POS],
        nB = dm[nB_POS],
        nit = dm[itr_POS],
        nbn = dm[bun_POS], 
        nth = dm[thn_POS],
        nS = dm[kp_POS],
        nR = dm[rpt_POS],
        tn = dm[tnit_POS],
        ntn = dm[ntn_POS];
    int i, j, k;
    double acc_pct[]={0,0,0}, *init, **sims,
        tnw = REAL(getListElement(x,"tune.weight"))[0];
    SEXP inits = getListElement(x,"inits"), ans, ans_tmp;

    //allocate memory for struct and simulated values
    bcpglm_str *da = (bcpglm_str *) R_alloc(1,sizeof(bcpglm_str)) ;
    da_parm *dap = (da_parm *) R_alloc(1,sizeof(da_parm)) ;
    da->dap = dap ;
    da->simTvec = ivect(nP) ;

    // fill in struct dap 
    dap->nO = nO ;
    dap->nP = nP;    
    dap->nB = nB ;
    dap->ygt0 = YPO_ELT(x) ;
    dap->Y= Y_ELT(x) ;
    dap->offset= OFFSET_ELT(x) ;
    dap->weights= PWT_ELT(x) ;
    dap->X = X_ELT(x) ;
    dap->link_power = LKP_ELT(x)[0] ;
    dap->beta = BETA_ELT(x) ;
    dap->phi = PHI_ELT(x)[0] ;
    dap->p = P_ELT(x)[0];
  
    // fill in struct da
    da->mu = MU_ELT(x) ;
    da->eta = ETA_ELT(x) ;
    da->pbeta_mean = REAL(getListElement(x,"pbeta.mean")) ;
    da->pbeta_var = REAL(getListElement(x,"pbeta.var")) ;
    da->mh_beta_var = REAL(getListElement(x,"ebeta.var"));
    da->mh_p_var = REAL(getListElement(x,"ep.var"))[0];
    da->mh_phi_var = REAL(getListElement(x,"ephi.var"))[0];
    da->bound_p = BDP_ELT(x) ;
    da->bound_phi = REAL(getListElement(x,"bound.phi"))[0] ;
    
    // update eta and mu
    cpglm_fitted(da->eta, da->mu, (double*) NULL, dap) ;

    // tune the scale parameter for M-H update    
    if (tn){
        if (nR>0)
            Rprintf("Tuning phase...\n");
        int etn = ceil(tn *1.0/ntn) ;  // # iters per tuning loop
        sims = dmatrix(etn,nB+2) ;
        double *beta_sims = dvect(etn*nB) ;
        double m_p =0,  sam_p_var=0;
        double m_phi =0,  sam_phi_var=0;
        double *sam_beta_var = dvect(nB*nB) ;

        for (k=0;k<ntn;k++) {
            bcpglm_mcmc_lat(da,  0, etn, 0, 1, sims, acc_pct);
            // adjust proposal variance for p and phi
            for (i=0;i<etn;i++){   // compute sample mean
                m_p += sims[i][nB+1] ;
                m_phi += sims[i][nB] ;
            }
            m_p /= etn ;
            m_phi /= etn ;
            for (i=0;i<etn;i++){   // compute sample covariance 
                sam_p_var += (sims[i][nB+1]-m_p)*(sims[i][nB+1]-m_p);
                sam_phi_var += (sims[i][nB]-m_phi)*(sims[i][nB]-m_phi);
            }
            sam_p_var /= etn-1 ;
            sam_phi_var /= etn-1 ;
            // adjust proposal covariance if accept percentage is too low or too high
            if (acc_pct[0]<0.4 || acc_pct[0] > 0.6)
                da->mh_p_var = tnw * da->mh_p_var + (1-tnw) * sam_p_var  ;
            if (acc_pct[2]<0.4 || acc_pct[2] > 0.6)
                da->mh_phi_var = tnw * da->mh_phi_var + (1-tnw) * sam_phi_var  ;

            // adjust vcov for beta            
            // compute sample vcov matrix
            for (i=0;i<etn;i++){
                for (j=0;j<nB;j++)
                    beta_sims[i+j*etn] = sims[i][j] ; // convert to long vector
            }
            cplm_cov(etn, nB, beta_sims, sam_beta_var) ;
            // adjust scale 
            if (acc_pct[1]<0.15 || acc_pct[1] > 0.35){
                for (i=0;i<nB*nB;i++)
                    da->mh_beta_var[i] = tnw * da->mh_beta_var[i] + (1-tnw) * sam_beta_var[i];
            }
        }
        if (nR>0){
            Rprintf("Acceptance rate in the last tuning phase:  beta(%4.2f%%), phi(%4.2f%%), p(%4.2f%%)\n",
                    acc_pct[1]*100, acc_pct[2]*100, acc_pct[0]*100);
            Rprintf("-----------------------------------------\n");
        }
     }
    // run Markov chains
    PROTECT(ans=allocVector(VECSXP,dm[chn_POS])) ;
    if (nR>0){
        Rprintf("Markov Chain Monte Carlo starts...\n");
        Rprintf("-----------------------------------------\n");
    }
    // simulations
    sims = dmatrix(nS,nB+2) ;
    for (k=0;k<dm[chn_POS];k++){
        if (nR>0)
            Rprintf("Start Markov chain %d\n", k+1);
        // re-initialize 
        init = REAL(VECTOR_ELT(inits,k));
        dap->beta = init ;
        dap->phi = init[nB] ;
        dap->p = init[nB+1];
        // update eta and mu
        cpglm_fitted(da->eta, da->mu, (double*) NULL, dap) ;
        bcpglm_mcmc_lat(da, nR, nit, nbn, nth, sims,acc_pct );
        //return result    
        PROTECT(ans_tmp=allocMatrix(REALSXP, nS, nB+2));
        for (j=0;j<nB+2;j++){		
            for (i=0;i<nS;i++)		
                REAL(ans_tmp)[i+nS*j]= sims[i][j] ;		
        }
        SET_VECTOR_ELT(ans, k, ans_tmp);
        UNPROTECT(1) ;
        if (nR>0)
            Rprintf("-----------------------------------------\n");
    }
    UNPROTECT(1) ;
    if (nR>0)
        Rprintf("Markov Chain Monte Carlo ends!\n");
    return ans ;
    
}

