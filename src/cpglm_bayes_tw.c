/************************************************************/
/*   Function for the Markov Chain Monte Carlo algorithm    */
/*    in the Compound Poisson Generalized Linear Model      */
/*        using direct tweedie density evaluations          */
/*              Author:  Wayne Zhang                        */
/*            actuary_zhang@hotmail.com                     */
/************************************************************/

/**
 * @file cpglm_bayes_tw.c
 * @brief Function for implementing the MCMC algorithm
 * in the Compound Poisson Generalized Linear Model using
 * direct tweedie density evaluation
 * @author Wayne Zhang                         
 */

#include "cplm.h"

/************************************************/
/*   Function to compute full conditionals      */  
/************************************************/

/**
 * posterior log density of the index parameter p
 *
 * @param x value of p at which the log density is to be calculated
 * @param data a void struct to work with arms, cocerced to SEXP internally
 *
 * @return log posterior density
 */
static double bcpglm_post_p_tw(double x, void *data){
    SEXP da= data;
    int *dm = DIMS_ELT(da) ;
    double *Y = Y_ELT(da), *mu = MU_ELT(da), phi = PHI_ELT(da)[0] ;
    return -0.5*dl2tweedie(dm[nO_POS], Y, mu, phi, x) ;
}


/**
 * posterior log density of the index parameter phi
 *
 * @param x value of phi at which the log density is to be calculated
 * @param data a void struct to work with arms, cocerced to SEXP internally
 *
 * @return log posterior density
 */
static double bcpglm_post_phi_tw(double x, void *data){
    SEXP da = data ;
    int *dm = DIMS_ELT(da) ;
    double *Y = Y_ELT(da), *mu = MU_ELT(da), p = P_ELT(da)[0] ;
    return -0.5*dl2tweedie(dm[nO_POS], Y, mu, x, p)  ;
}

/**
 * posterior log density of of the vector of beta
 *
 * @param x vector of values for beta
 * @param data void struct that is coerced to SEXP
 *
 * @return log posterior density for beta
 */
static double bcpglm_post_beta_tw(double *x,  void *data){
    SEXP da = data ;
    int *dm = DIMS_ELT(da) ;
    int nO = dm[nO_POS],
        nP = dm[nP_POS],
        nB = dm[nB_POS];    
    int i, kk, *ygt0 = YPO_ELT(da) ;
    double ld=0, p= P_ELT(da)[0], phi = PHI_ELT(da)[0];
    double p2=2-p, p1=p-1;
    double *offset= OFFSET_ELT(da), *wts =PWT_ELT(da), *X = X_ELT(da),
        *Y = Y_ELT(da),*link_power = LKP_ELT(da),
        *eta = ETA_ELT(da), *mu = MU_ELT(da);
    double *pbeta_mean = REAL(getListElement(da,"pbeta.mean")),
        *pbeta_var = REAL(getListElement(da,"pbeta.var")) ;

    // update mu
    cplm_eta(eta, nO, nB, X, x, (double *) NULL, offset);
    cplm_mu_eta(mu, (double *) NULL, nO, eta, *link_power) ;
    
    // loglikelihood from data
    for (i=0; i<nO; i++)
        ld += pow(mu[i],p2) * wts[i];
    ld /= (- phi*p2) ;
    for (i=0; i<nP; i++){
        kk = ygt0[i] ;
        ld += - Y[kk]*pow(mu[kk],-p1)*wts[kk] /(phi*p1);
    }
    // prior info
    for (i=0;i<nB;i++)
        ld += -0.5*(x[i]-pbeta_mean[i])*(x[i]-pbeta_mean[i])/pbeta_var[i] ;
    return ld ;
}

/************************************************/
/*     Main function to fit compound Poisson    */
/*     GLM using Monte Carlo Markov Chains      */
/************************************************/
/**
 * implement MCMC for compound Poisson GLM using tweedie density evaluation
 *
 * @param da a list object
 * @param nR report interval
 * @param nit number iterations
 * @param nbn number of burn-in
 * @param nth thinning rate
 * @param sims a 2d array to store simulation results
 * @param acc_pct a vector of length 3 to store acceptance percentage
 *
 */

static void bcpglm_mcmc_tw(SEXP da, int nR, int nit, int nbn, int nth,
                           double **sims, double *acc_pct){
    int *dm = DIMS_ELT(da) ;
    int nO = dm[nO_POS],
        nB = dm[nB_POS];
    int i, j, iter,  ns ;
    int  acc=0, accept[]={0,0,0};
    // bound for p and phi
    double xl_p = BDP_ELT(da)[0], xr_p =BDP_ELT(da)[1],
        xr_phi=REAL(getListElement(da,"bound.phi"))[0];
    // proposal covariance matrix
    double *mh_beta_var = REAL(getListElement(da,"ebeta.var")),
        mh_p_var = REAL(getListElement(da,"ep.var"))[0],
        mh_phi_var = REAL(getListElement(da,"ephi.var"))[0];
    double *offset= OFFSET_ELT(da),*X = X_ELT(da),
        *link_power = LKP_ELT(da), *eta = ETA_ELT(da), *mu = MU_ELT(da),
        *beta= BETA_ELT(da), *p = P_ELT(da), *phi= PHI_ELT(da) ;
    double xtemp, *beta_sim=Calloc(nB, double) ;
    double p_sd = sqrt(mh_p_var), phi_sd = sqrt(mh_phi_var) ;

    // update eta and mu
    cplm_eta(eta, nO, nB, X, beta, (double *) NULL, offset);
    cplm_mu_eta(mu, (double *) NULL, nO, eta, *link_power) ;

    GetRNGstate() ;
    for (iter=0;iter<nit;iter++){
        if (nR>0 && (iter+1)%nR==0)
            Rprintf("Iteration: %d \n ", iter+1) ;
        R_CheckUserInterrupt() ;
        
        // M-H update of p using truncated normal
        acc = metrop_tnorm_rw(*p, p_sd, xl_p, xr_p, &xtemp, 
                                bcpglm_post_p_tw, (void *) da);	
        *p = xtemp ;
        accept[0] += acc ;
        R_CheckUserInterrupt() ;
        
        //Metropolis-Hasting block update              
        acc = metrop_mvnorm_rw(nB, beta, mh_beta_var,
                               beta_sim, bcpglm_post_beta_tw, (void *)da) ;
        Memcpy(beta, beta_sim, nB) ;
        accept[1] += acc ;
    
        // update eta and mu
        cplm_eta(eta, nO, nB, X, beta, (double *) NULL, offset);
        cplm_mu_eta(mu, (double *) NULL, nO, eta, *link_power) ;
        R_CheckUserInterrupt() ;
        
        // M-H update of phi using truncated normal
        acc = metrop_tnorm_rw(*phi, phi_sd, 0, xr_phi, &xtemp, 
                                bcpglm_post_phi_tw, (void *) da);
        *phi = xtemp ;
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
                sims[ns][j] = beta[j];
            sims[ns][nB] = *phi  ;
            sims[ns][nB+1] = *p ;      
            } 
    }
    PutRNGstate() ;
    Free(beta_sim);
    // compute acceptance percentage
    for (i=0;i<3;i++)
        acc_pct[i] = accept[i]*1.0/nit ;
}

/**
 * implement MCMC for compound Poisson GLM
 *
 * @param da a list object
 *
 * @return the simulated values
 *
 */

SEXP bcpglm_gibbs_tw (SEXP da){
    // get dimensions
    int *dm = DIMS_ELT(da) ;
    int nB = dm[nB_POS],
        nit = dm[itr_POS],
        nbn = dm[bun_POS], 
        nth = dm[thn_POS],
        nS = dm[kp_POS],
        nR = dm[rpt_POS],
        tn = dm[tnit_POS],
        ntn = dm[ntn_POS];
    int i, j, k;
    double acc_pct[]={0,0,0}, *init, **sims,
        tnw = REAL(getListElement(da,"tune.weight"))[0];
    SEXP inits = getListElement(da,"inits"), ans, ans_tmp;
    
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
        // proposal covariance matrix 
        double *mh_beta_var = REAL(getListElement(da,"ebeta.var")),
            *mh_p_var = REAL(getListElement(da,"ep.var")),
            *mh_phi_var = REAL(getListElement(da,"ephi.var"));

        for (k=0;k<ntn;k++) {
            bcpglm_mcmc_tw(da,  0, etn, 0, 1, sims, acc_pct);
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
                *mh_p_var = tnw * (*mh_p_var) + (1-tnw) * sam_p_var  ;
            if (acc_pct[2]<0.4 || acc_pct[2] > 0.6)
                *mh_phi_var = tnw * (*mh_phi_var) + (1-tnw) * sam_phi_var  ;

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
                    mh_beta_var[i] = tnw * mh_beta_var[i] + (1-tnw) * sam_beta_var[i];
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
        Memcpy(BETA_ELT(da),init, nB) ;
        PHI_ELT(da)[0] = init[nB] ;
        P_ELT(da)[0] = init[nB+1];
        bcpglm_mcmc_tw(da, nR, nit, nbn, nth, sims, acc_pct);
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
