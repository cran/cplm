#include "utilities.h"

/* This is the struct needed for computing the hessian matrix
   nY: number of positive values
   nO: number of observations
   nS: number of simulated samples for T_i
   nB: length of regression parameters
   Y: vector of POSITIVE response  	       
   X: design matrix for beta
   simT: array of simulated T's
   et: vector of expected values of T's
   w: vector of weights
   sw: sum of weights (w)
   ygt0: vector of index of positive values - used for mu
   linkpower: power of link function
*/

typedef struct {
  int nY, nO, nS, nB, k, *ygt0;	
  double *Y, **X, **simT, *et, *w,  sw, linkpower, *offset ;
}  hessGlm;

// eta^{-1} to compute mu
double linkInv(double x, double linkpower){
    double xx ;
    if (linkpower==0)
        xx = exp(x);
    else
        xx= pow(x,1/linkpower) ;
    return xx ;           
}

// d mu / d eta to compute gradient
double muEta(double x, double linkpower){
    double xx ;
    if (linkpower==0)
        xx = exp(x);
    else
        xx= (1/linkpower)*pow(x, 1/linkpower-1) ;
    return xx ;           
}

/*
  function to compute the mean (mu) of cpglm
*/
void muGlm(double *theta, void *ex, double *mu,
             double *eta, int etaT){
  hessGlm *da= ex ;
  double xb ;
  int i, j ;

  for (i=0; i< da->nO; i++) {
        xb =0 ;
        for (j=0;j<da->nB;j++){
            xb += da->X[i][j] * theta[j] ;
        }
	xb += da->offset[i] ;
        if (etaT==1)
            eta[i] = xb ;        
        mu[i] = linkInv(xb,da->linkpower);
    }
}

/************************************************/
/*  Function to compute variance of gradient    */
/*   For each t_i, compute gradient, then       */
/*   compute expectation and second moments     */
/************************************************/

/*
  function to compute the log posterior of joint density
  given the k_{th} simulated T_i
*/

double dJLogLikGlm1(double *theta, void *ex){
    hessGlm *da = ex ;
    int i,  k=da->k ;
    double lp=0, phi=theta[da->nB], p=theta[da->nB+1];
    double p2=2-p, p1=p-1;
    double *mu = vect(da->nO) ;
    // compute mu 
    muGlm(theta, da, mu, NULL, 0) ;

    // compute log likelihood
    for (i=0;i<da->nO;i++) 
        lp +=pow(mu[i],p2) ;    
    lp *= - (1/(phi*p2)) ;
    for (i=0; i<da->nY; i++){ 
        lp += - da->Y[i]* pow(mu[da->ygt0[i]],-p1)/(phi*p1)
             - lgammafn(da->simT[i][k]*p2/p1)
             + da->simT[i][k]*(p2/p1*(log(da->Y[i])-log(p1))-log(phi)/p1 - log(p2));
    }
    return lp ;
}

/*
  compute  derivative:
   - beta and phi is computed analytically
   - p is computed numerically 
*/ 

void dJLogLikGlmGrad1 (int n, double *par, double *gr,  void *ex){
    hessGlm *da =ex ;
    double *mu = vect(da->nO) ;
    double *eta = vect(da->nO) ;
    double *mu1p= vect(da->nO) ;
    double *mueta =vect(da->nO) ;    
    int i, j ;
    double phi=par[da->nB], p=par[da->nB+1];
    double p2=2-p, p1=p-1;

    // compute mu 
    muGlm(par, da, mu, eta, 1) ;

    // pre-compute mu^(1-p) and d mu/ d eta for later use
    for (i=0;i<da->nO;i++){
        mu1p[i]= pow(mu[i],-p1) ;
        mueta[i] = muEta(eta[i],da->linkpower) ;
    }

    // derivative for beta
    for (j=0; j<n-2;j++){
        gr[j] = 0 ;
        for (i=0;i<da->nO;i++)
            gr[j] += - mu1p[i] * mueta[i] * da->X[i][j];
        for (i=0;i<da->nY;i++)
            gr[j] += da->Y[i]* mu1p[da->ygt0[i]]/mu[da->ygt0[i]]
                * mueta[da->ygt0[i]] * da->X[da->ygt0[i]][j] ;
        gr[j] *= (1/phi) ;
    }

    // derivative for phi
    gr[n-2]=0 ;
    for (i=0;i<da->nO;i++)
        gr[n-2] += mu1p[i]*mu[i];
    gr[n-2] *= 1/ (phi*phi*p2) ;
    for (i=0;i<da->nY;i++) {
        gr[n-2] += da->Y[i]*mu1p[da->ygt0[i]]/(p1*phi*phi)
            - da->simT[i][da->k]/(p1*phi) ;
    }

    // derivative for p, using numerical solution
    double dp1, dp2 ;
    par[n-1] += EPS ;
    dp1 = dJLogLikGlm1(par,da) ;
    par[n-1] -= 2*EPS ;
    dp2 = dJLogLikGlm1(par,da) ;
    gr[n-1]=(dp1-dp2)/(2*EPS);
    par[n-1] += EPS ;        
}

// compute hessian numerically
void dJLogLikGlmHess1 (int n, double *par, double **hess,  void *ex){
    hessGlm *da =ex ;
    int i, j ;    
    double *df1 = vect(n);
    double *df2 = vect(n);
    
    for (i = 0; i < n; i++) {
	par[i] +=  EPS;
	dJLogLikGlmGrad1(n, par, df1, da);
	par[i] -=  2 * EPS;
        dJLogLikGlmGrad1(n, par, df2, da);
	for (j = 0; j < n; j++)
	    hess[i][j] = (df1[j] - df2[j])/ (2*EPS) ;
	par[i] +=EPS;
    }
}

/************************************************/
/*  Function to compute second derivative       */
/*   of E(L), i.e. the hessian of Q             */
/************************************************/

/*
  function to compute the expectation of the
  log posterior of joint density
*/

double dJLogLikGlm2(double *theta, void *ex){
    hessGlm *da = ex ;
    int i, j ;
    double lp=0, phi=theta[da->nB], p=theta[da->nB+1];
    double p2=2-p, p1=p-1,  elgt;
    double *mu = vect(da->nO) ;

    // compute mu 
    muGlm(theta, da, mu, NULL, 0) ;

    // compute expected Logliklihood
    for (i=0; i<da->nO; i++)
        lp += pow(mu[i],p2);
    lp *= - (1/(phi*p2)) ;      
    for (i=0; i<da->nY; i++){
         // compute E(lgamma(t*p2/p1))
        elgt =0.0 ; 
        for (j=0; j<da->nS;j++){
            // weighted by w
            elgt += lgammafn(da->simT[i][j]*p2/p1)* da->w[j] ;  
        }
        elgt *= 1/da->sw ;
        lp += - da->Y[i]*pow(mu[da->ygt0[i]],-p1)/(phi*p1)
                - elgt
            + da->et[i]*(p2/p1*(log(da->Y[i])-log(p1))-log(phi)/p1 - log(p2));
    }    
    return lp ;
}

/*
  compute  derivative:
   - beta and phi is computed analytically
   - p is computed numerically 
*/ 

void dJLogLikGlmGrad2 (int n, double *par, double *gr,  void *ex){
    hessGlm *da =ex ;
    double *mu = vect(da->nO) ;
    double *eta = vect(da->nO) ;
    double *mu1p= vect(da->nO) ;
    double *mueta =vect(da->nO) ;    
    int i, j ;
    double phi=par[da->nB], p=par[da->nB+1];
    double p2=2-p, p1=p-1;
    
    // compute mu 
    muGlm(par, da, mu, eta, 1) ;

    // pre-compute mu^(1-p) and d mu/ d eta for later use
    for (i=0;i<da->nO;i++){
        mu1p[i]= pow(mu[i],-p1) ;
        mueta[i] = muEta(eta[i],da->linkpower) ;
    }

    // derivative for beta
    for (j=0; j<n-2;j++){
        gr[j] = 0 ;
        for (i=0;i<da->nO;i++)
            gr[j] += - mu1p[i] * mueta[i] * da->X[i][j];
        for (i=0;i<da->nY;i++)
            gr[j] += da->Y[i]* mu1p[da->ygt0[i]]/mu[da->ygt0[i]]
                * mueta[da->ygt0[i]] * da->X[da->ygt0[i]][j] ;
        gr[j] *= (1/phi) ;
    }

    // derivative for phi
    gr[n-2]=0 ;
    for (i=0;i<da->nO;i++)
        gr[n-2] += mu1p[i]*mu[i];
    gr[n-2] *= 1/ (phi*phi*p2) ;
    for (i=0;i<da->nY;i++) {
        gr[n-2] += da->Y[i]*mu1p[da->ygt0[i]]/(p1*phi*phi)
            - da->et[i]/(p1*phi) ;
    }

    // derivative for p, using numerical solution
    double dp1, dp2 ;
    par[n-1] += EPS ;
    dp1 = dJLogLikGlm2(par,da) ;
    par[n-1] -= 2*EPS ;
    dp2 = dJLogLikGlm2(par,da) ;
    gr[n-1]=(dp1-dp2)/(2*EPS);
    par[n-1] += EPS ;    
}

// compute hessian numerically
void dJLogLikGlmHess2 (int n, double *par, double **hess,  void *ex){
    hessGlm *da =ex ;
    int i, j ;    
    double *df1 = vect(n);
    double *df2 = vect(n);

    for (i = 0; i < n; i++) {
	par[i] +=  EPS;
	dJLogLikGlmGrad2(n, par, df1, da);
	par[i] -=  2 * EPS;
        dJLogLikGlmGrad2(n, par, df2, da);
	for (j = 0; j < n; j++)
	    hess[i][j] = (df1[j] - df2[j])/ (2*EPS) ;
	par[i] +=EPS;
    }
}

/*
  Main function:
  * compute hessian using Monte Carlo methods
  - X, design matrix
  - Y, the vector of positive responses
  - ygt0, vector of index for postive y's
  - theta, all the parameters in model: beta, phi,p in order
  - simT, matrix of simulated T's
  - w, vector of weights
  - linkpower, power for link function
  - n, number of observations
  * return value is a vector of length 2
*/

SEXP hessGlmEst (SEXP X, SEXP Y, SEXP ygt0,
              SEXP theta, SEXP simT, SEXP w,
		 SEXP linkpower, SEXP n, SEXP offset){    
    // dimensions of return matrix
    int nY = LENGTH(Y),
      nS=LENGTH(w),
        nO= INTEGER(n)[0],
        nB = LENGTH(theta) - 2 ;
    int i, j, k;
    SEXP ans, ans1, ans2, ans3, names ;
    hessGlm *da ;

    // allocate memory 
    da = (hessGlm *) R_alloc(1,sizeof(hessGlm)) ;
    da->Y = vect(nY) ;
    da->ygt0 =(int *) R_alloc(nY,sizeof(int)) ;
    da->et = vect(nY) ;
    da->w = vect(nS) ;
    da->offset=vect(nO) ;
    da->simT = matrix(nY,nS) ;
    da->X = matrix(nO,nB) ;

    // fill in struct da
    da->nY = nY;
    da->nO = nO ;
    da->nS = nS ;
    da->nB = nB ;
    da->k =0 ;

    da->linkpower= REAL(linkpower)[0] ;
    for (i=0;i<nY;i++){
        da->Y[i]=REAL(Y)[i] ;
        // ygt0 from R does not start from 0
        da->ygt0[i]=INTEGER(ygt0)[i]-1 ;
    }
    for (j=0;j<nS;j++){
        da->w[j] = REAL(w)[j] ;
        for (i=0;i<nY;i++){
            da->simT[i][j]=INTEGER(simT)[i+nY*j] ;
        }
    }
    for (j=0;j<nB;j++){
        for (i=0;i<nO;i++){
            da->X[i][j]=REAL(X)[i+nO*j] ;
        }
    }    
    for (i=0;i<nO;i++)
      da->offset[i]=REAL(offset)[i] ;
    da->sw = cumsum(da->w,nS) ;    
    // compute weighted E(T_i)
    for (i=0;i<nY;i++){
        da->et[i]=cumwsum(da->simT[i],da->w,nS)/da->sw;
    } 

    /************************************************/
    /* compute var of gradient= E(L*L') - E(L)*E(L)' */
    /************************************************/
    // E(L*L')
    double **VL = matrix(nB+2,nB+2) ;
    double **ELLt = matrix(nB+2,nB+2) ;
    double **ELELt = matrix(nB+2,nB+2) ;
    double *EL = vect(nB+2) ;
    double *grad = vect(nB+2) ;

    //initialize EL and ELLt
    for (j=0;j<nB+2;j++){
        EL[j] = 0 ;
        for (k=0;k<nB+2;k++){
            ELLt[j][k] = 0 ;
        }
    }

    // compute E(L) and E(L*L') 
    for (i=0;i<nS;i++){
        da ->k=i ;
        dJLogLikGlmGrad1(nB+2,REAL(theta), grad, da);
        for (j=0;j<nB+2;j++){
            // compute E(L)
            EL[j] += grad[j]/nS ;
            for (k=0;k<nB+2;k++){
                // compute E(L*L')
                ELLt[j][k] += grad[j]*grad[k]/nS;
            }
        }
    }

    // compute E(L)*E(L')
    for (i=0;i<nB+2;i++){
        for (j=0;j<nB+2;j++){
            ELELt[i][j]=EL[i]*EL[j];
            VL[i][j]=ELLt[i][j]-ELELt[i][j] ;
        }
    }

    /************************************************/
    /*          compute hessian of Q:               */
    /************************************************/
    double **hQ = matrix(nB+2, nB+2) ;
    dJLogLikGlmHess2 (nB+2, REAL(theta), hQ, da);    
        
    /************************************************/
    /*        return the above two as a list        */
    /************************************************/

    PROTECT(ans=allocVector(VECSXP,3)) ;
    // VL 
    PROTECT(ans1=allocMatrix(REALSXP,nB+2,nB+2));
    // dQ
    PROTECT(ans2=allocMatrix(REALSXP,nB+2,nB+2));
    // VQ
    PROTECT(ans3=allocMatrix(REALSXP,nB+2,nB+2));

    for (j=0;j<nB+2;j++){
        for (i=0;i<nB+2;i++){
            REAL(ans1)[i+(nB+2)*j]=VL[i][j] ;
            REAL(ans2)[i+(nB+2)*j]=hQ[i][j] ;
	    // variance of Q, used for computing normal 
	    // confidence interval to increase sample size
            REAL(ans3)[i+(nB+2)*j]=ELLt[i][j]/nS ;
        }
    }
    SET_VECTOR_ELT(ans, 0, ans1);
    SET_VECTOR_ELT(ans, 1, ans2);
    SET_VECTOR_ELT(ans, 2, ans3);
    // set names attributes
    PROTECT(names = allocVector(STRSXP, 3));
    SET_STRING_ELT(names, 0, mkChar("VL"));
    SET_STRING_ELT(names, 1, mkChar("HQ"));
    SET_STRING_ELT(names, 2, mkChar("VQ"));
    setAttrib(ans, R_NamesSymbol, names);
    UNPROTECT(5) ;
    return ans ;

}





