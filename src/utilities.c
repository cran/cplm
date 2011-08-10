//  
//  Some utility functions used in cplm
//       Wayne Zhang
//

#include "utilities.h"

// allocate memory for a double vector of size n
double * vect(int n)
{
    return (double *)R_alloc(n, sizeof(double));
}

// allocate memory for a double matrix of nr \times nc
double ** matrix(int nr, int nc)
{
    int i;
    double **m;
    m = (double **) R_alloc((nr + 1), sizeof(double *));
    for (i = 0; i <= nr; i++)
	m[i] = (double*) R_alloc((nc + 1), sizeof(double));
    return m;
}

// cumulative sum of a vector
double cumsum(double *x, int n){
    int i;
    double sm=0 ;
    for (i=0;i<n;i++)
        sm +=x[i] ;
    return sm ;
}

// weight cumulative sum of a vector by weight w
double cumwsum(double *x, double *w, int n){
    int i;
    double sm=0 ;
    for (i=0;i<n;i++)
        sm +=x[i]*w[i] ;
    return sm ;
}

/*
 wrapper for the univariate lbfgsb function:
 - x, the value of the parameter
 - lower, the lower bound
 - upper, the upper bound
 - val, the value of the function
 - fn, function to be minimized
 - gr, derivative of function
 - ex, struct to be passed to fn and gr
*/

void lbfgsbU(double *x, double lower, double upper, double *val,
             optimfn fn, optimgr gr, void *ex,  int *conv){
    // nbd needed for lbfgsb 
    int nbd ;
    if (!R_FINITE(lower)) {
        if (!R_FINITE(upper)) nbd = 0; else nbd = 3;
    } else {
        if (!R_FINITE(upper)) nbd = 1; else nbd = 2;
    }

    // default the parameters needed for lbfgsb 
    double factr=1e7, pgtol=0 ;
    int fncount, grcount;
    int lmm=5, maxit=1000, trace=0, nREPORT=10 ;
    char msg[60] ;

    lbfgsb(1,lmm,x, &lower, &upper, &nbd, val, fn, gr,
           conv, (void *)ex, factr, pgtol, &fncount, &grcount,
           maxit,msg ,trace,nREPORT) ;
} 


// extract element from a list
SEXP getListElement (SEXP list, char *str)
{
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;
  for (i = 0; i < length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  return elmt;
}
