/************************************************************/
/*   Function for the Markov Chain Monte Carlo algorithm    */
/*    in the Compound Poisson Generalized Linear Mixed      */
/*    Model using direct tweedie density evaluations        */
/*              Author:  Wayne Zhang                        */
/*            actuary_zhang@hotmail.com                     */
/************************************************************/

/**
 * @file bcpglmm_tw.c
 * @brief Function for implementing the MCMC algorithm
 * in the Compound Poisson Generalized Linear Mixed Model using
 * direct tweedie density evaluation
 * @author Wayne Zhang                         
 */

#include "cplm.h"
#include "Matrix.h"		 /* for cholmod functions */

/** cholmod_common struct initialized  */
extern cholmod_common c;

/*
 * initiate cholmod and set the factorization form to be LL'
 */
/*
SEXP init(){
M_R_cholmod_start(&c);
    c.final_ll = 1;
    return R_NilValue ;
    }
SEXP finish(){
    M_cholmod_finish(&c);
    return R_NilValue ;
    }
*/

/************************************************/
/*             Some utility functions           */  
/************************************************/

/**
 * Compute the mean in cpglmm
 *
 * @param da a list object
 *
 */
static void cpglmm_fitted(SEXP da){    
    int *dm = DIMS_ELT(da) ;
    int nO = dm[nO_POS], nB = dm[nB_POS], i1 = 1  ;
    double *offset= OFFSET_ELT(da), *X = X_ELT(da),
        *link_power = LKP_ELT(da), *beta = BETA_ELT(da),
        *eta = ETA_ELT(da), *mu = MU_ELT(da), one[] = {1,0};
    CHM_DN ceta, u = AS_CHM_DN(getListElement(da,"u"));
    CHM_SP Zt = Zt_ELT(da);
    R_CheckStack();
    // update eta
    Memcpy(eta, offset, nO);
    // eta := eta + X * beta 
    F77_CALL(dgemv)("N", &nO, &nB, one, X, &nO,
		    beta, &i1, one, eta, &i1);
    ceta = N_AS_CHM_DN(eta, nO, 1);
    R_CheckStack();
    if (!M_cholmod_sdmult(Zt, 1 , one, one, u, ceta, &c))
        error(_("cholmod_sdmult error returned"));
    // update mu
    cplm_mu_eta(mu, (double *) NULL, nO, eta, *link_power) ;
}

/**
 * Set parameter to the k_th initial values provided in the inits slot
 *
 * @param da a list object
 * @param k indicates the k_th set of initial values
 *
 */
static void set_init(SEXP da, int k){
    int *dm = DIMS_ELT(da) ;
    int i, pos = 0, nB = dm[nB_POS], nU = dm[nU_POS],
        nT = dm[nT_POS], *nc = NCOL_ELT(da);
    SEXP inits = getListElement(da, "inits"),
        Sig = getListElement(da, "Sigma");
    double *Sigi, *init = REAL(VECTOR_ELT(inits,k));
    Memcpy(BETA_ELT(da),init, nB) ;
    PHI_ELT(da)[0] = init[nB] ;
    P_ELT(da)[0] = init[nB+1];
    Memcpy(U_ELT(da),init+nB+2, nU) ;
    for (i=0;i<nT;i++){
        Sigi = REAL(VECTOR_ELT(Sig,i)) ;
        Memcpy(Sigi, init+nB+2+nU+pos, nc[i]*nc[i]) ;
        pos += nc[i]*nc[i] ;
    }    
}

/**
 * Set parameter to the ns_th column of the simulation results
 *
 * @param da a list object
 * @param ns indicates the ns_th column
 * @param sims matrix to store simulations results 
 *
 */
static void set_sims(SEXP da, int ns, double **sims){
    SEXP Sig = getListElement(da, "Sigma");
    int *dm = DIMS_ELT(da) ;
    int i, pos = 0, nB = dm[nB_POS], nU = dm[nU_POS],
        nT = dm[nT_POS], *nc = NCOL_ELT(da);
    double *Sigi ;
    for (i=0;i<nB;i++)
        sims[ns][i] = BETA_ELT(da)[i];
    sims[ns][nB] = PHI_ELT(da)[0] ;
    sims[ns][nB+1] = P_ELT(da)[0] ;
    for (i=0;i<nU;i++)
        sims[ns][nB+2+i] = U_ELT(da)[i];
    for (i=0;i<nT;i++){
        Sigi = REAL(VECTOR_ELT(Sig, i));
        Memcpy(sims[ns]+nB+2+nU+pos,Sigi, nc[i]*nc[i]) ;
        pos += nc[i]*nc[i] ;
        } 
}

/************************************************/
/*   Function to compute full conditionals      */  
/************************************************/

/**
 * posterior log density of the index parameter p
 *
 * @param x value of p at which the log density is to be calculated
 * @param data a void struct, cocerced to SEXP internally
 *
 * @return log posterior density
 */

double bcpglmm_post_p_tw(double x, void *data){
    SEXP da= data;
    int *dm = DIMS_ELT(da) ;
    double *Y = Y_ELT(da), *mu = MU_ELT(da), phi = PHI_ELT(da)[0] ;
    return -0.5*dl2tweedie(dm[nO_POS], Y, mu, phi, x) ;
}

/**
 * posterior log density of the index parameter phi
 *
 * @param x value of phi at which the log density is to be calculated
 * @param data a void struct, cocerced to SEXP internally
 *
 * @return log posterior density
 */

static double bcpglmm_post_phi_tw(double x, void *data){
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
static double bcpglmm_post_beta_tw(double *x,  void *data){
    SEXP da = data ;
    int *dm = DIMS_ELT(da) ;
    int nO = dm[nO_POS],
        nP = dm[nP_POS],
        nB = dm[nB_POS];    
    int i, kk, *ygt0 = YPO_ELT(da) ;
    double ld=0, p= P_ELT(da)[0], phi = PHI_ELT(da)[0];
    double p2=2-p, p1=p-1;
    double *wts =PWT_ELT(da), *Y = Y_ELT(da), *mu = MU_ELT(da), 
        *pbeta_mean = PBM_ELT(da), *pbeta_var = PBV_ELT(da),
        *beta= BETA_ELT(da), *beta_old = Alloca(nB, double) ;
    R_CheckStack() ;
    
    // update mu
    Memcpy(beta_old, beta, nB) ;
    Memcpy(beta, x, nB) ;
    cpglmm_fitted(da) ;
    Memcpy(beta, beta_old, nB) ;
    
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

/**
 * posterior log density of of the vector of u
 *
 * @param x vector of values for beta
 * @param data void struct that is coerced to SEXP
 *
 * @return log posterior density for u
 */
static double bcpglmm_post_u_tw(double *x,  void *data){
    SEXP da = data,
        Sig = getListElement(da,"Sigma") ;
    int *dm = DIMS_ELT(da) ;
    int nO = dm[nO_POS], nP = dm[nP_POS],
        nU = dm[nU_POS], nT = dm[nT_POS];       
    int i, j, k, kk,
        *ygt0 = YPO_ELT(da), *Gp = Gp_ELT(da),
        *nc = NCOL_ELT(da), *nlev= NLEV_ELT(da);
    int mc = imax(nc, nT); 
    double ld=0, p= P_ELT(da)[0], phi = PHI_ELT(da)[0];
    double p2=2-p, p1=p-1, *Sigi, *wts =PWT_ELT(da),
        *Y = Y_ELT(da), *u=U_ELT(da), *mu = MU_ELT(da),
        *xv = Alloca(mc, double), *iv = Alloca(mc*mc, double),
        *u_old = Alloca(nU, double);
    R_CheckStack() ;
    
    // update mu
    Memcpy(u_old, u, nU) ;
    Memcpy(u, x, nU);
    cpglmm_fitted(da) ;
    Memcpy(u, u_old, nU) ;
    
    // loglikelihood from data
    for (i=0; i<nO; i++)
        ld += pow(mu[i],p2) * wts[i];
    ld /= (- phi*p2) ;
    for (i=0; i<nP; i++){
        kk = ygt0[i] ;
        ld += - Y[kk]*pow(mu[kk],-p1)*wts[kk] /(phi*p1);
    }
    
    // prior info    
    for (i=0;i<nT;i++){
        Sigi = REAL(VECTOR_ELT(Sig, i));        
        solve_po(nc[i],Sigi, iv) ;
        for (j=0;j<nlev[i];j++){
            for (k=0;k<nc[i];k++){
                kk = Gp[i]+ j + k*nlev[i] ;
                xv[k]= x[kk] ; // u vector when nc >1
            }
            ld += dmvnorm(nc[i], xv, (double*) NULL, iv);
        }            
    }
    return ld ;
}


/************************************************/
/*     Main function to fit compound Poisson    */
/*     GLMM using Monte Carlo Markov Chains     */
/************************************************/
/**
 * implement MCMC for compound Poisson GLMM using tweedie density evaluation
 *
 * @param da a list object
 * @param nR report interval
 * @param nit number iterations
 * @param nbn number of burn-in
 * @param nth thinning rate
 * @param sims a 2d array to store simulation results
 * @param acc_pct a vector of length 4 to store acceptance percentage
 *
 */
  
static void bcpglmm_mcmc_tw(SEXP da, int nR, int nit, int nbn, int nth,
                           double **sims, double *acc_pct){
    SEXP Sig = getListElement(da,"Sigma"),
        pSig = getListElement(da, "pSigma"), pSigi ;
    int *dm = DIMS_ELT(da) ;
    int nB = dm[nB_POS], nT = dm[nT_POS], nU = dm[nU_POS];
    int i, j, iter,  ns, pos, acc=0, *Gp = Gp_ELT(da),
        *nc = NCOL_ELT(da), *nlev= NLEV_ELT(da), accept[]={0,0,0,0};
    int mc = imax(nc, nT);

    // bound for p and phi
    double xl_p = BDP_ELT(da)[0], xr_p =BDP_ELT(da)[1],
        xr_phi = BDPHI_ELT(da)[0];
    // proposal covariance matrix etc..
    double *mh_beta_var = EBV_ELT(da), *mh_u_var = EUV_ELT(da),
        mh_p_var = EPV_ELT(da)[0], mh_phi_var = EPHIV_ELT(da)[0],        
        *beta= BETA_ELT(da), *u = U_ELT(da), *p = P_ELT(da),
        *phi= PHI_ELT(da),
        *beta_sim = Alloca(nB, double), *beta_m = Alloca(nB, double),
        *u_sim=Alloca(nU, double), *u_m = Alloca(nU, double),
        *scl = Alloca(mc*mc, double), *scl2 = Alloca(mc*mc, double);
    double xtemp, su, *Sigi, p_sd = sqrt(mh_p_var), phi_sd = sqrt(mh_phi_var);
    R_CheckStack() ;
    
    // update eta and mu
    cpglmm_fitted(da) ;
    
    GetRNGstate() ;
    for (iter=0;iter<nit;iter++){
        if (nR>0 && (iter+1)%nR==0)
            Rprintf("Iteration: %d \n ", iter+1) ;
        R_CheckUserInterrupt() ;
       
        // M-H update of p using truncated normal
        acc = metrop_tnorm_rw(*p, p_sd, xl_p, xr_p, &xtemp, 
                                bcpglmm_post_p_tw, (void *) da);	
        *p = xtemp ;
        accept[0] += acc ;
        R_CheckUserInterrupt() ;
          
        //Metropolis-Hasting block update of beta        
        Memcpy(beta_m, beta, nB);  // need this because of the copying step in post_beta               
        acc = metrop_mvnorm_rw(nB, beta_m, mh_beta_var,
                               beta_sim, bcpglmm_post_beta_tw, (void *)da) ;                
        Memcpy(beta, beta_sim, nB) ;
        accept[1] += acc ;
        cpglmm_fitted(da);        
        R_CheckUserInterrupt() ;
        
        //Metropolis-Hasting block update of u
        Memcpy(u_m, u, nU);     
        acc = metrop_mvnorm_rw(nU, u_m, mh_u_var,
                               u_sim, bcpglmm_post_u_tw, (void *)da) ;
        Memcpy(u, u_sim, nU) ;
        accept[2] += acc ;
        cpglmm_fitted(da) ;        
        R_CheckUserInterrupt() ;
        
        // M-H update of phi using truncated normal
        acc = metrop_tnorm_rw(*phi, phi_sd, 0, xr_phi, &xtemp, 
                                bcpglmm_post_phi_tw, (void *) da);
        *phi = xtemp ;
        accept[3] += acc ;
        R_CheckUserInterrupt() ;
        
        // direct simulation of Sigma due to conjugacy
        for (i=0;i<nT;i++){
            Sigi = REAL(VECTOR_ELT(Sig, i));
            pSigi = VECTOR_ELT(pSig, i) ;
            if (nc[i]==1){
                // simulate from bounded inverse-Gamma
                su = norm(u+Gp[i], nlev[i]) ;                    
                Sigi[0] = 1 / rgamma(0.5*nlev[i]+IGSHP_ELT(pSigi)[0],
                                     1/(su*su*0.5+IGSCL_ELT(pSigi)[0])) ;                
            }
            else {                   
                // simulate from inverse-Wishart
                mult_xtx(nlev[i],nc[i], u+Gp[i], scl) ;  // t(x) * (x)
                pos = 0;
                for (j=0;j<nc[i]*nc[i];j++)
                    scl[j] += IWSCL_ELT(pSigi)[j] ;
                solve_po(nc[i], scl, scl2) ;
                rwishart(nc[i], (double) nlev[i]+ IWDF_ELT(pSigi)[0], scl2, scl) ;
                solve_po(nc[i], scl, Sigi) ;                  
            }
        }
        
        // print out acceptance rate if necessary
        if (nR>0 && (iter+1)%nR==0){
            Rprintf(_("Acceptance rate: beta(%4.2f%%), u(%4.2f%%), phi(%4.2f%%), p(%4.2f%%),\n"),
                    accept[1]*1.0/(iter+1)*100, accept[2]*1.0/(iter+1)*100,
                    accept[3]*1.0/(iter+1)*100, accept[0]*1.0/(iter+1)*100 );
                    } 
        // store results        
        if (iter>=nbn &&  (iter+1-nbn)%nth==0 ){
            ns = (iter+1-nbn)/nth -1;
            set_sims(da, ns, sims) ;
        }
    }
    PutRNGstate() ;
    // compute acceptance percentage
    for (i=0;i<4;i++)
        acc_pct[i] = accept[i]*1.0/nit ;
}

/**
 * tune the proposal covariance matrix
 *
 * @param da an input list object
 * @param acc_pct a vector to store the acceptance rate
 *
 */
static void bcpglmm_tune_tw(SEXP da, double *acc_pct){
    int *dm = DIMS_ELT(da) ;
    int nB = dm[nB_POS], nA = dm[nA_POS],
        nU = dm[nU_POS], nR = dm[rpt_POS],
        tn = dm[tnit_POS], ntn = dm[ntn_POS];
    int i, j, k;
    double **sims, tnw = REAL(getListElement(da,"tune.weight"))[0];

    if (nR>0)
        Rprintf("Tuning phase...\n");
    int etn = ceil(tn *1.0/ntn) ;  // # iters per tuning loop
    sims = dmatrix(etn,nA) ;
    // store samples and sample covariance
    double *p_sims = dvect(etn), *phi_sims = dvect(etn),
        *beta_sims = dvect(etn*nB), *u_sims = dvect(etn*nU);
    double sam_p_var,  sam_phi_var, 
        *sam_beta_var = dvect(nB*nB), *sam_u_var = dvect(nU*nU);
    // proposal covariance matrix 
    double *mh_beta_var = EBV_ELT(da), *mh_u_var = EUV_ELT(da),
        *mh_p_var = EPV_ELT(da), *mh_phi_var = EPHIV_ELT(da);
    
    for (k=0;k<ntn;k++) {
        // run mcmc
        bcpglmm_mcmc_tw(da,  0, etn, 0, 1, sims, acc_pct);                
        // convert to long vector
        for (i=0;i<etn;i++){   
            p_sims[i] = sims[i][nB+1] ;
            phi_sims[i] = sims[i][nB] ;
            for (j=0;j<nB;j++)
                beta_sims[i+j*etn] = sims[i][j] ; 
            for (j=0;j<nU;j++)
                u_sims[i+j*etn] = sims[i][nB+2+j] ;                

        }
        // adjust proposal variance for p and phi
        cov(etn,1,p_sims, &sam_p_var) ;
        cov(etn,1,phi_sims, &sam_phi_var) ;        
        if (acc_pct[0]<0.4 || acc_pct[0] > 0.6)
            *mh_p_var = tnw * (*mh_p_var) + (1-tnw) * sam_p_var  ;
        if (acc_pct[3]<0.4 || acc_pct[3] > 0.6)
            *mh_phi_var = tnw * (*mh_phi_var) + (1-tnw) * sam_phi_var  ;

        // adjust vcov for beta and u
        cov(etn, nB, beta_sims, sam_beta_var) ;
        cov(etn, nU, u_sims, sam_u_var) ;
        if (acc_pct[1]<0.15 || acc_pct[1] > 0.35){
            for (i=0;i<nB*nB;i++)
                mh_beta_var[i] = tnw * mh_beta_var[i] + (1-tnw) * sam_beta_var[i];
        }
        if (acc_pct[2]<0.15 || acc_pct[2] > 0.35){
            for (i=0;i<nU*nU;i++)
                mh_u_var[i] = tnw * mh_u_var[i] + (1-tnw) * sam_u_var[i];
        }        
    }
    if (nR>0){
        Rprintf("Acceptance rate in the last tuning phase:  beta(%4.2f%%), u(%4.2f%%), phi(%4.2f%%), p(%4.2f%%)\n",
                acc_pct[1]*100, acc_pct[2]*100, acc_pct[3]*100, acc_pct[0]*100);
        Rprintf("-----------------------------------------\n");
    }
}

/**
 * implement MCMC for compound Poisson GLMM using direct density evaluation
 *
 * @param da an input list object
 *
 * @return the simulated values
 *
 */

SEXP bcpglmm_gibbs_tw (SEXP da){
    // get dimensions
    int *dm = DIMS_ELT(da) ;
    int nA = dm[nA_POS],
        nit = dm[itr_POS], nbn = dm[bun_POS], 
        nth = dm[thn_POS], nS = dm[kp_POS],
        nR = dm[rpt_POS], tn = dm[tnit_POS];
    int i, j, k;
    double acc_pct[]={0,0,0,0},**sims;
    SEXP ans, ans_tmp;
    
    // tune the scale parameter for M-H update    
    if (tn)
        bcpglmm_tune_tw(da, acc_pct) ;
        
    // run Markov chains
    PROTECT(ans=allocVector(VECSXP,dm[chn_POS])) ;
    if (nR>0){
        Rprintf("Markov Chain Monte Carlo starts...\n");
        Rprintf("-----------------------------------------\n");
    }
    // simulations
    sims = dmatrix(nS,nA) ;
    for (k=0;k<dm[chn_POS];k++){
        if (nR>0)
            Rprintf("Start Markov chain %d\n", k+1);
        // re-initialize
        set_init(da, k) ;
        bcpglmm_mcmc_tw(da, nR, nit, nbn, nth, sims, acc_pct);
        //return result    
        PROTECT(ans_tmp=allocMatrix(REALSXP, nS, nA));
        for (j=0;j<nA;j++){		
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
    