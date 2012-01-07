/************************************************************/
/*   Function for the Markov Chain Monte Carlo algorithm    */
/*    in the Compound Poisson Generalized Linear Model      */
/*            using the latent variable approach            */
/*              Author:  Wayne Zhang                        */
/*            actuary_zhang@hotmail.com                     */
/************************************************************/

/**
 * @file bcpglm_lat.c
 * @brief Function for implementing the MCMC algorithm
 * in the Compound Poisson Generalized Linear Model using
 * the latent variable approach
 * @author Wayne Zhang                         
 */

#include "cplm.h"


/************************************************/
/*   Function to compute full conditionals      */  
/************************************************/

/**
 * posterior log density of the index parameter p
 * (this function is the same as that in bcpglmm_lat)
 *
 * @param x value of p at which the log density is to be calculated
 * @param data a void struct, cocerced to SEXP internally
 *
 * @return log posterior density
 */
static double bcpglm_post_p_lat(double x, void *data){
    SEXP da = data ;
    double p_old = P_ELT(da)[0] ;
    P_ELT(da)[0] = x ;
    double lp = cplm_llik_lat(da) ;
    P_ELT(da)[0] = p_old ;
    return lp ;
}

/**
 * posterior log density of the dispersion parameter phi
 * (this function is the same as that in bcpglmm_lat)
 *
 * @param x value of phi at which the log density is to be calculated
 * @param data a void struct, cocerced to SEXP internally
 *
 * @return log posterior density
 */
static double bcpglm_post_phi_lat(double x, void *data){
    SEXP da = data ;
    int *dm = DIMS_ELT(da), *ygt0 = YPO_ELT(da),
        *simT = SIMT_ELT(da)  ;
    int i, kk, nO = dm[nO_POS], nP = dm[nP_POS] ;
    double *Y = Y_ELT(da), *mu = MU_ELT(da), *wts =PWT_ELT(da),
        p = P_ELT(da)[0];
    double ld=0, p2=2-p, p1= p-1 ;
    for (i=0; i<nO; i++)
        ld += pow(mu[i],p2) * wts[i];
    ld /= (- x*p2) ;
    for (i=0; i<nP; i++){
        kk = ygt0[i] ;
        ld += - Y[kk]*pow(mu[kk],-p1)*wts[kk] /(x*p1) 
            - log(x)* simT[i]/p1;
    }
    return ld ;
}


/**
 * posterior log density of of the vector of beta
 *
 * @param x vector of values for beta
 * @param data void struct that is coerced to SEXP internally
 *
 * @return log posterior density for beta
 */
double bcpglm_post_beta_lat(double *x,  void *data){
    SEXP da = data ;
    int *dm = DIMS_ELT(da) ;
    int nO = dm[nO_POS], nP = dm[nP_POS], nB = dm[nB_POS];    
    int i, kk, *ygt0 = YPO_ELT(da) ;
    double ld=0, p= P_ELT(da)[0], phi = PHI_ELT(da)[0];
    double p2=2-p, p1=p-1;
    double *wts =PWT_ELT(da), *Y = Y_ELT(da), *mu = MU_ELT(da), 
        *pbeta_mean = PBM_ELT(da), *pbeta_var = PBV_ELT(da) ;
    
    // update mu
    cpglm_fitted_x(x, da) ;
    
    // loglikelihood from data
    for (i=0; i<nO; i++)
        ld += pow(mu[i],p2) * wts[i];
    ld /= (- phi*p2) ;
    for (i=0; i<nP; i++){
        kk = ygt0[i] ;
        ld += - Y[kk]* pow(mu[kk],-p1)* wts[kk] /(phi*p1);
    }
    // prior info
    for (i=0;i<nB;i++)
        ld += -0.5*(x[i]-pbeta_mean[i])* (x[i]-pbeta_mean[i])/pbeta_var[i] ;
    return ld ;
}


/************************************************/
/*     Main function to fit compound Poisson    */
/*     GLM using Monte Carlo Markov Chains      */
/************************************************/

/**
 * MCMC simulation for compound Poisson GLM
 *
 * @param da a SEXP struct
 * @param nR report interval
 * @param nit number of iterations
 * @param nbn number of burn-ins
 * @param nth thinning rate
 * @param sims a 2d array to store simulation results
 * @param acc_pct acceptance percentage 
 *
 */
static void bcpglm_mcmc_lat(SEXP da, int nR, int nit, int nbn,
                            int nth, double **sims, double *acc_pct){
    int *dm = DIMS_ELT(da), *kk = K_ELT(da) ;
    int nB = dm[nB_POS], nP=dm[nP_POS];
    int  acc, accept[]={0,0,0}; 
    int i, j, iter, ns ;

    // proposal covariance matrix etc..
    double *mh_beta_var = EBV_ELT(da), mh_p_var = EPV_ELT(da)[0],
        mh_phi_var = EPHIV_ELT(da)[0], *beta= BETA_ELT(da),
        *p = P_ELT(da), *phi= PHI_ELT(da),
        *beta_sim = Alloca(nB, double) ;
    R_CheckStack() ;
    double xtemp, xl_p = BDP_ELT(da)[0], xr_p =BDP_ELT(da)[1],
        xr_phi = BDPHI_ELT(da)[0], p_sd = sqrt(mh_p_var),
        phi_sd= sqrt(mh_phi_var) ;
    
    GetRNGstate() ;
    for (iter=0;iter<nit;iter++){
        if (nR>0 && (iter+1)%nR==0)
            Rprintf("Iteration: %d \n ", iter+1) ;
        R_CheckUserInterrupt() ;
        
        // update latent variable T using rejection sampling
        for (i=0;i<nP;i++){    
            *kk = i ;
            cplm_rlatT_reject(da) ;
        }
        R_CheckUserInterrupt() ;
        
        // M-H update of p using truncated normal
        acc = metrop_tnorm_rw(*p, p_sd, xl_p, xr_p, &xtemp, 
                              bcpglm_post_p_lat, (void *) da);
        *p = xtemp ;
        accept[0] += acc ;
        R_CheckUserInterrupt() ;
        
        //Metropolis-Hasting block update              
        acc = metrop_mvnorm_rw(nB, beta, mh_beta_var,
                               beta_sim, bcpglm_post_beta_lat, (void *)da) ;
        Memcpy(beta, beta_sim, nB) ;
        accept[1] += acc ;
        cpglm_fitted(da) ;
        R_CheckUserInterrupt() ;
    

        // M-H update of phi using truncated normal
        acc = metrop_tnorm_rw(*phi, phi_sd, 0, xr_phi, &xtemp, 
                              bcpglm_post_phi_lat, (void *) da);
        *phi = xtemp ;
        accept[2] += acc ;
        R_CheckUserInterrupt() ;
        
        // print out acceptance rate if necessary
        if (nR>0 && (iter+1)%nR==0){
            Rprintf(_("Acceptance rate: beta(%4.2f%%), phi(%4.2f%%), p(%4.2f%%)\n"),
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
    // compute acceptance percentage
    for (i=0;i<3;i++)
        acc_pct[i] = accept[i]*1.0/nit ;
}


/**
 * implement MCMC for compound Poisson GLM using latent variables
 *
 * @param da a list object
 *
 * @return the simulated values 
 *
 */

SEXP bcpglm_gibbs_lat (SEXP da){
    // get dimensions
    int *dm = DIMS_ELT(da) ;
    int nB = dm[nB_POS], nit = dm[itr_POS],
        nbn = dm[bun_POS], nth = dm[thn_POS],
        nS = dm[kp_POS], nR = dm[rpt_POS],
        tn = dm[tnit_POS], ntn = dm[ntn_POS];
    int i, j, k;
    double acc_pct[]={0,0,0}, *init, **sims,
             tnw = REAL(getListElement(da,"tune.weight"))[0];
    SEXP inits = getListElement(da,"inits"), ans, ans_tmp;
      
    // update eta and mu
    cpglm_fitted(da) ;

    // tune the scale parameter for M-H update    
    if (tn){
        int etn = ceil(tn *1.0/ntn) ;  // # iters per tuning loop
        double *beta_sims = dvect(etn*nB), *p_sims = dvect(etn),
            *phi_sims = dvect(etn), sam_p_var, sam_phi_var,
            *sam_beta_var = dvect(nB*nB) ;
        double *mh_beta_var = EBV_ELT(da), *mh_p_var = EPV_ELT(da),
            *mh_phi_var = EPHIV_ELT(da) ;
        sims = dmatrix(etn,nB+2) ;
    
        if (nR>0)
            Rprintf("Tuning phase...\n");
    
        for (k=0;k<ntn;k++) {
            bcpglm_mcmc_lat(da,  0, etn, 0, 1, sims, acc_pct);
            // convert to long vector
            for (i=0;i<etn;i++){   
                p_sims[i] = sims[i][nB+1] ;
                phi_sims[i] = sims[i][nB] ;
                for (j=0;j<nB;j++)
                    beta_sims[i+j*etn] = sims[i][j] ; 
            }
            // adjust proposal variance for p and phi
            cov(etn,1,p_sims, &sam_p_var) ;
            cov(etn,1,phi_sims, &sam_phi_var) ;        
            if (acc_pct[0]<0.4 || acc_pct[0] > 0.6)
                *mh_p_var = tnw * (*mh_p_var) + (1-tnw) * sam_p_var  ;
            if (acc_pct[2]<0.4 || acc_pct[2] > 0.6)
                *mh_phi_var = tnw * (*mh_phi_var) + (1-tnw) * sam_phi_var  ;

            // adjust vcov for beta
            cov(etn, nB, beta_sims, sam_beta_var) ;
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
        // update eta and mu
        cpglm_fitted(da) ;
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

