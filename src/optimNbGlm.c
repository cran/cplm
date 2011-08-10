#include "utilities.h"

/* This is the struct needed for optimization
   nY: number of positive values
   nO: number of observations
   nS: number of simulated samples for T_i
   Y: vector of POSITIVE response  	       
   mu: vector of means
   simT: array of simulated T's
   et: vector of expected values of T's
   w: vector of weights
   phi: scale parameter
   p: index parameter for Tweedie
   sw: sum of weights (w)
   ygt0: vector of index of positive values - used for mu
*/

typedef struct {
  int nY, nO, nS ;
  double phi,p, sw;	
  double *Y,  *mu, **simT, *et, *w ;
  int *ygt0;  
}  emCpglm;

/*
  function to compute the log posterior of p
*/

double dLogPostP(double x, void *ex){
    emCpglm *da = ex ;
    int i, j ;
    double lp=0, p2=2-x, p1=x-1, elgt ;

    for (i=0; i<da->nO; i++)
        lp += pow(da->mu[i],p2);
    lp *= - (1/(da->phi*p2)) ;
    for (i=0; i<da->nY; i++){
         // compute E(lgamma(t*p2/p1))
        elgt =0.0 ; 
        for (j=0; j<da->nS;j++){
            elgt += lgammafn(da->simT[i][j]*p2/p1)* da->w[j] ;  // weighted by w        
        }
        elgt *= 1/da->sw ;
        lp += - da->Y[i]*pow(da->mu[da->ygt0[i]],-p1)/(da->phi*p1)
                - elgt
            + da->et[i]*(p2/p1*(log(da->Y[i])-log(p1))-log(da->phi)/p1 - log(p2));
    }
    return lp ;

}

/*
  The following two functions are to compute the log posterior
  of p and its gradient, as required in lfgsbU
*/

static double dLogPostPFunc(int n, double *par, void *ex){
    emCpglm *da =ex ;    
    double ans= dLogPostP(*par, da) * SCALE;
    return ans ;
}

// compute numerical derivative
static void dLogPostPGrad (int n, double *par, double *gr,  void *ex){
    emCpglm *da =ex ;
    double par1 = *par +EPS, par2 =*par-EPS ;
    *gr=(dLogPostP(par1,da)-dLogPostP(par2,da))/(2*EPS) * SCALE;
}

/*
  Main function:
  *  optimization for phi and p
  - mu, vector of means for all obs
  - Y, the vector of positive responses
  - ygt0, vector of index for postive y's
  - phi, scale parameter
  - p, index parameter
  - simT, matrix of simulated T's
  - w, vector of weights
  - bd, lower and upper bound for p
  * return value is a vector of length 2
*/
SEXP optimNbGlm ( SEXP Y, SEXP mu, SEXP ygt0,
              SEXP p, SEXP phi, SEXP simT, SEXP w,
             SEXP bd){    
    // dimensions of return matrix
    int nY = LENGTH(Y), nO=LENGTH(mu), nS=LENGTH(w) ;    
    double pnew, dpval,phinew;
    int i, j, conv;
    SEXP ans ;
    emCpglm *da ;

    // allocate memory for struct and Y, mu, simT
    da = (emCpglm *) R_alloc(1,sizeof(emCpglm)) ;
    da->Y = vect(nY) ;
    da->mu = vect(nO) ;
    da->ygt0 =(int *) R_alloc(nY,sizeof(int)) ;
    da->et = vect(nY) ;
    da->w = vect(nS) ;
    da->simT = matrix(nY,nS) ;
    
    // fill in struct da
    da->nY = nY;
    da->nO = nO ;
    da->nS = nS ;
    for (i=0;i<nY;i++){
        da->Y[i]=REAL(Y)[i] ;
        // ygt0 from R does not start from 0
        da->ygt0[i]=INTEGER(ygt0)[i]-1 ;
    }
    for (i=0;i<nO;i++)
        da->mu[i]=REAL(mu)[i] ;
    for (j=0;j<nS;j++){
        da->w[j] = REAL(w)[j] ;
        for (i=0;i<nY;i++){
            da->simT[i][j]=INTEGER(simT)[i+nY*j] ;
        }
    }
    da->phi = REAL(phi)[0] ;
    da->p = REAL(p)[0];
    da->sw = cumsum(da->w,nS) ;
    // compute weighted E(T_i)
    for (i=0;i<nY;i++){
        da->et[i]=cumwsum(da->simT[i],da->w,nS)/da->sw;
    } 

    // optimization and return results
    PROTECT(ans=allocVector(REALSXP,2));

    // optimization of phi- direct computation
    phinew=0 ;
    double p1=REAL(p)[0]-1, p2=2-REAL(p)[0];
    for (i=0;i<nO;i++)
        phinew += pow(da->mu[i],p2) ;
    phinew *= p1/p2 ;
    for (i=0;i<nY;i++)
        phinew += da->Y[i]*pow(da->mu[da->ygt0[i]],-p1);
    phinew *= 1/cumsum(da->et,nY) ;
    REAL(ans)[0] = phinew ;
    // update phi
    da->phi = phinew ;  
  
    // optimization of p
    pnew = REAL(p)[0] ;
    lbfgsbU(&pnew, REAL(bd)[0], REAL(bd)[1], &dpval,
            dLogPostPFunc, dLogPostPGrad, (void *)da,  &conv) ;

    // if converged, set to new value; o/w, retain
    REAL(ans)[1] = (conv > 0) ? REAL(p)[0] : pnew ;
    UNPROTECT(1) ;
    return ans ;
}
