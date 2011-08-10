#include "utilities.h"

/* This is the struct needed for optimization 
   Y: vector of POSITIVE response  	       
   k: the index for the current T or Y
   phi: scale parameter
   p: index parameter for Tweedie
   lambda: original mean for the truncated poisson proposal
*/

typedef struct {
  int k ;
  double phi, p, lambda;	
  double *Y ;	
}  gpOptim;

/* 
   Log posterior of the k_th T
   - x, a point at which the density to be computed
   - y, the positive response value corresponding to the t_i
   - phi, scale parameter
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
double gpDiffFunc (int n, double *par,  void *data){
  gpOptim *da = data ;
  double ldiff = dLogPostT(par[0],da->Y[da->k],da->phi, da->p)
                   - dLogTruncPois(par[0],da->lambda) ;
  return ldiff ;
}

/*
  derivative of the gpDiffFunc 
*/
void gpDiffGrad (int n, double *par, double *gr,  void *data){
    gpOptim *da =data ;
    double p1 = da->p -1 ;
    double p2 = 2 - da->p ;
    *gr = p2/p1*log(da->Y[da->k]/p1)-log(da->phi)/p1-log(p2)-
        log(da->lambda)-p2/p1*digamma(p2/p1*par[0]);
}

// the following two functions scale the above two by -1 for max
double gpDiffFunc2 (int n, double *par,  void *data){
  double ldiff = gpDiffFunc(n, par, data) * SCALE ;
  return ldiff ;
}

void gpDiffGrad2 (int n, double *par, double *gr,  void *data){
  gpDiffGrad(n,par,gr,data) ;
  *gr = *gr * SCALE ;
}

/*
  Rejections sampling of the latent variable T  *
  - n, number of samples need for each T_i
  - Y, the vector of positive responses
  - phi, scale parameter
  - p, index parameter
  * proposal is the zero truncated poisson with
     orignal mean mathed to the approximate mode
     of the posterior of T
  * return value is a length(Y) \times n matrix of
     simulated values
*/
SEXP sampleTRej (SEXP n, SEXP Y, SEXP phi, SEXP p){
    // dimensions of return matrix
    int nc = INTEGER(n)[0] ;
    int nr = LENGTH(Y) ;
    double par, val, pb, u, xtemp2;
    int i, j, xtemp, accept, conv ;
    SEXP ans ;
    gpOptim *da ;

    // allocate memory for struct and Y
    da = (gpOptim *) R_alloc(1,sizeof(gpOptim)) ;
    da->Y = vect(nr) ;

    // fill in struct da
    for (i=0;i<nr;i++)
        da->Y[i]=REAL(Y)[i] ;
    da->phi = REAL(phi)[0] ;
    da->p = REAL(p)[0] ;

    PROTECT(ans=allocMatrix(INTSXP,nr,nc));
    GetRNGstate() ;
    for (i=0; i<nr;i++){
        // initial value of par
        par = 2 ;
        da->k=i ;
        da->lambda = lambdaTruncPois(da->Y[da->k],da->phi, da->p) ;
        // find t that results in the largest difference
        // between target and proposal
        lbfgsbU(&par, 1, INFINITY, &val,gpDiffFunc2, gpDiffGrad2,
                (void *)da,  &conv) ;
	val *= 1/SCALE ;
        for (j=0;j<nc;j++){
            accept =0 ;
            while (accept==0) {
                xtemp =rTruncPois(da->lambda);
                // gpDiffFunc only takes double!
                xtemp2 = xtemp ;
                pb = exp(gpDiffFunc(1, &xtemp2,da)-val);
                u = runif(0,1);
                if (u<=pb){
                    INTEGER(ans)[i+nr*j] = xtemp ;
                    accept =1 ;
                }
            }
        }
    }
    PutRNGstate() ;
    UNPROTECT(1) ;
    return ans ;
}


/*
  Compute weights in importance sampling 
 */
SEXP importWeight(SEXP nc, SEXP Y, SEXP sT, SEXP parm){
  int nS = INTEGER(nc)[0], nY = LENGTH(Y) ;
  int i, j ;
  double **simT;
  SEXP parm1, parm2, ans ;
  parm1 = getListElement(parm, "old");
  parm2 = getListElement(parm, "new");
  simT = matrix(nY,nS) ;

  // assign T_i
   for (j=0;j<nS;j++){
        for (i=0;i<nY;i++){
           simT[i][j]=INTEGER(sT)[i+nY*j] ;
        }
    }
   // compute weight
   PROTECT(ans=allocVector(REALSXP,nS)) ;
   for (j=0;j<nS;j++){
     REAL(ans)[j] =0.0 ;
     for (i=0;i<nY;i++){
       REAL(ans)[j] += dLogPostT(simT[i][j],REAL(Y)[i], REAL(parm2)[0], REAL(parm2)[1])
	 -dLogPostT(simT[i][j],REAL(Y)[i], REAL(parm1)[0], REAL(parm1)[1]);
     }
     REAL(ans)[j] = exp(REAL(ans)[j]) ;
   }
   UNPROTECT(1) ;
   return ans ;
}


