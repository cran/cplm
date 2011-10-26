/************************************************************/
/*   Header files for C functions in the cplm package       */
/*              Author:  Wayne Zhang                        */
/*            actuary_zhang@hotmail.com                     */
/************************************************************/

/**
 * @file cplm.h
 * @brief header files to be included in the cplm C functions
 * @author Wayne Zhang                         
*/


#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <R_ext/Applic.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>


/** constant used in numerical gradient and hessian */
#define EPS (0.001)
/** constant used in optim to convert min to max */
#define SCALE (-1)   

/** zero an array */
#define AZERO(x, n) {int _I_, _SZ_ = (n); for(_I_ = 0; _I_ < _SZ_; _I_++) (x)[_I_] = 0;}


#ifndef DA_PARM_STRUCT 
#define DA_PARM_STRUCT
/** struct to store data and parameters */ 
typedef struct {
    int nO ;              /**< number of Observations */
    int nP ;              /**< # positive observations */
    int nB ;              /**< # model cofficients*/
    int k ;               /**< index used internally by various functions */
    int *ygt0 ;           /**< row index of positive values */
    double *X ;           /**< design matrix */
    double *Y ;           /**< reponse variable */ 
    double *offset ;      /**< vector of offsets */
    double *weights ;     /**< prior weights */
    double *beta ;        /**< model coefficients */
    double phi ;          /**< dispersion parameter */
    double p ;            /**< index parameter */
    double link_power ;   /**< power of link function, as in tweedie */
    double lambda ;	  /**< original mean for the truncated Poisson proposal    */
}  da_parm ;

#endif

#ifdef ENABLE_NLS		/** Allow for translation of error messages */
#include <libintl.h>
#define _(String) dgettext ("cplm", String)
#else
#define _(String) (String)
#endif



/** positions in the dims vector in Bayesian models */
enum dimB {
    nO_POS=0,			/**<number of observations */
    nB_POS,			/**<number of coefficients */
    nP_POS,			/**<number of positive obs */
    chn_POS,			/**<number of chains */
    itr_POS,			/**<number of iterations */
    bun_POS,			/**<number of burn-in */
    thn_POS,			/**<number of thinning */
    kp_POS,			/**<number of simulations to keep in each chain */
    sims_POS,			/**<total number of simulations */
    rpt_POS,			/**<report frequency */
    tnit_POS,			/**<number of iterations used for tuning */
    ntn_POS			/**<number of tuning loops */    
};


// R list 
SEXP getListElement (SEXP list, char *str);

/**
 * Extract the element named nm from the list object obj and return a null pointer
 * if the element has length 0 or a pointer to the REAL contents.
 *
 * @param obj pointer to a list object
 * @param str pointer to a symbol naming the element to extract
 *
 * @return pointer to the REAL contents, if nonzero length, otherwise
 * a NULL pointer
 *
 */
static R_INLINE double *ELT_REAL_NULL(SEXP obj, char *str)
{
    SEXP pt = getListElement(obj, str);
    return LENGTH(pt) ? REAL(pt) : (double*) NULL; 
}

/** Return the double pointer to the X slot */
#define X_ELT(x) ELT_REAL_NULL(x, "X")

/** Return the double pointer to the y slot */
#define Y_ELT(x) ELT_REAL_NULL(x, "y")

/** Return the double pointer to the offset slot or (double*) NULL if
 * offset has length 0) */
#define OFFSET_ELT(x) ELT_REAL_NULL(x, "offset")

/** Return the double pointer to the pWt slot or (double*) NULL if
 * pWt has length 0) */
#define PWT_ELT(x) ELT_REAL_NULL(x, "pWt")

/** Return the integer pointer to the ygt0 slot or (double*) NULL if
 * pWt has length 0) */
#define YPO_ELT(x) INTEGER(getListElement(x, "ygt0"))

/** Return the integer pointer to the dims slot */
#define DIMS_ELT(x) INTEGER(getListElement(x, "dims"))

/** Return the double pointer to the fixef slot */
#define BETA_ELT(x) ELT_REAL_NULL(x, "beta")

/** Return the double pointer to the eta slot */
#define ETA_ELT(x) ELT_REAL_NULL(x, "eta")

/** Return the double pointer to the mu slot */
#define MU_ELT(x) ELT_REAL_NULL(x, "mu")

/** Return the double pointer to the muEta slot or (double*) NULL if
 * muEta has length 0) */
#define MUETA_ELT(x) ELT_REAL_NULL(x, "muEta")

/** Return the double pointer to the var slot or (double*) NULL if
 * var has length 0) */
#define VAR_ELT(x) ELT_REAL_NULL(x, "var")

/** Return the double pointer to the resid slot */
#define RESID_ELT(x) ELT_REAL_NULL(x, "resid")

/** Return the double pointer to the p slot or (double*) NULL if
 *  sqrtXWt has length 0) */
#define P_ELT(x) ELT_REAL_NULL(x, "p")

/** Return the double pointer to the phi slot or (double*) NULL if
 *  sqrtXWt has length 0) */
#define PHI_ELT(x) ELT_REAL_NULL(x, "phi")

/** Return the double pointer to the link_power slot or (double*) NULL if
 *  sqrtXWt has length 0) */
#define LKP_ELT(x) ELT_REAL_NULL(x, "link.power")

/** Return the double pointer to the bound_p slot  */
#define BDP_ELT(x) ELT_REAL_NULL(x,"bound.p")



// memory allocation utilities
double * dvect(int n) ;
double ** dmatrix(int nr, int nc) ;
int * ivect(int n) ;
int ** imatrix(int nr, int nc) ;

// simple computation
double dcumsum(double *x, int n) ;
int icumsum(int *x, int n) ;
double dcumwsum(double *x, double *w, int n) ;
double icumwsum(int *x, double *w, int n) ;
double norm (double *x, int n);
double dist (double *x, double *y, int n);

// univariate optimization 
void lbfgsbU(double *x, double lower, double upper, double *val,
             optimfn fn, optimgr gr, void *ex,  int *conv) ;

// glm related computation
double varFun (double mu, double p);
double linkFun(double mu, double link_power);
double linkInv(double eta, double link_power);
double mu_eta(double eta, double link_power) ;

void cplm_eta(double *eta, int nO, int nB, double *X,
              double *beta, double *b, double *offset);       
void cplm_mu_eta(double* mu, double* muEta, int n, double* eta, 
		   double link_power);
void cplm_varFun(double* var, double* mu, int n, double p);
void cpglm_fitted(double *eta, double *mu, double *muEta,
                 da_parm *dap);


// cplm
double cplm_post_latT(double x, double y, double phi, double p) ; 
double cplm_lambda_tpois(double y, double phi, double p) ;
// function to simulate from 0 truncated poisson
int cplm_rtpois(double lambda ) ;
// function to compute the density of 0-truncated poisson on log scale
double cplm_dtpois(double x, double lambda) ;
void cplm_rlatT_reject (int nS, int *ans, da_parm *dap) ;

void cplm_tw_glm(double *x, double *y, double *off, double *wts, 
		 double *beta, double vp, double lp, int n, int p) ;

double cplm_llikS(double *mu, double phi, double p,
		      int *simT, da_parm *dap);


// M-H 
int metrop_mvnorm_rw(int d, double *m, double *v, double *sn, 
		     double (*myfunc)(double *x, void *data), 
		     void *data);

int metrop_tnorm_rw( double m, double sd, double lb, double rb, double *sn, 
		     double (*myfunc)(double x, void *data), 
		     void *data) ;
void cplm_cov(int n, int p, double *x, double *ans) ;


// tweedie 
double dtweedie(double y, double mu,
                 double phi, double p);
double dl2tweedie(int n, double *y, double *mu,
                  double phi, double p) ;
