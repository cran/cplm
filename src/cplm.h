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
#include "Matrix.h"

/** zero an array */
#define AZERO(x, n) {int _I_, _SZ_ = (n); for(_I_ = 0; _I_ < _SZ_; _I_++) (x)[_I_] = 0;}

/* When appropriate, alloca is cleaner than malloc/free.  The storage
 * is freed automatically on return from a function. When using gcc the
 * builtin version is much faster. */
#ifdef __GNUC__
# undef alloca
# define alloca(x) __builtin_alloca((x))
#else
/* this is necessary (and sufficient) for Solaris 10: */
# ifdef __sun
#  include <alloca.h>
# endif
#endif

/** alloca n elements of type t */
#define Alloca(n, t)   (t *) alloca( (size_t) ( (n) * sizeof(t) ) )


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
    nT_POS,                     /**<number of terms in cpglmm */
    nU_POS,                     /**<number of random coefficients */
    nA_POS,                     /**<number of parameters in the variance components */
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

/** Allocate (alloca) a cholmod_sparse struct, populate it with values
 * from the Zt slot and return the pointer. */
#define Zt_ELT(x) AS_CHM_SP(getListElement(x, "Zt"))

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

/** Return the double pointer to the ranef slot */
#define U_ELT(x) ELT_REAL_NULL(x, "u")

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

/** Return the double pointer to the link.power slot or (double*) NULL if
 *  sqrtXWt has length 0) */
#define LKP_ELT(x) ELT_REAL_NULL(x, "link.power")

/** Return the double pointer to the bound.p slot  */
#define BDP_ELT(x) ELT_REAL_NULL(x,"bound.p")

/** Return the double pointer to the bound.phi slot  */
#define BDPHI_ELT(x) ELT_REAL_NULL(x,"bound.phi")

/** Return the double pointer to the pbeta.mean slot  */
#define PBM_ELT(x) ELT_REAL_NULL(x,"pbeta.mean")

/** Return the double pointer to the pbeta.var slot  */
#define PBV_ELT(x) ELT_REAL_NULL(x,"pbeta.var")

/** Return the integer pointer to the Gp slot or (double*) NULL if
 * Gp has length 0) */
#define Gp_ELT(x) INTEGER(getListElement(x, "Gp"))

/** Return the double pointer to the ebeta.var slot  */
#define EBV_ELT(x) ELT_REAL_NULL(x,"ebeta.var")

/** Return the double pointer to the eu.var slot  */
#define EUV_ELT(x) ELT_REAL_NULL(x,"eu.var")

/** Return the double pointer to the ephi.var slot  */
#define EPHIV_ELT(x) ELT_REAL_NULL(x,"ephi.var")

/** Return the double pointer to the ep.var slot  */
#define EPV_ELT(x) ELT_REAL_NULL(x,"ep.var")

/** Return the integer pointer to the ncol slot  */
#define NCOL_ELT(x) INTEGER(getListElement(x, "ncol"))

/** Return the integer pointer to the nlev slot */
#define NLEV_ELT(x) INTEGER(getListElement(x, "nlev"))

/** Return the double pointer to the igamma.shape slot  */
#define IGSHP_ELT(x) ELT_REAL_NULL(x,"igamma.shape")

/** Return the double pointer to the igamma.scale slot  */
#define IGSCL_ELT(x) ELT_REAL_NULL(x,"igamma.scale")

/** Return the double pointer to the iwish.scale slot  */
#define IWSCL_ELT(x) ELT_REAL_NULL(x,"iwish.scale")

/** Return the double pointer to the iwish.df slot  */
#define IWDF_ELT(x) ELT_REAL_NULL(x,"iwish.df")

/** Return the integer pointer to the simT slot */
#define SIMT_ELT(x) INTEGER(getListElement(x, "simT"))

/** Return the integer pointer to the k slot */
#define K_ELT(x) INTEGER(getListElement(x, "k"))

/** Return the double pointer to the lambda slot  */
#define LAM_ELT(x) ELT_REAL_NULL(x,"lambda")


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
double dmax (double *x, int n) ;
int imax (int *x, int n) ;
    
// univariate optimization 
void lbfgsbU(double *x, double lower, double upper, double *val,
             optimfn fn, optimgr gr, void *ex,  int *conv) ;

// glm related computation
double varFun (double mu, double p);
double linkFun(double mu, double link_power);
double linkInv(double eta, double link_power);
double mu_eta(double eta, double link_power) ;

void cpglm_eta(double *eta, int nO, int nB, double *X,
              double *beta,  double *offset);       
void cplm_mu_eta(double* mu, double* muEta, int n, double* eta, 
		   double link_power);
void cplm_varFun(double* var, double* mu, int n, double p);
void cpglm_fitted(SEXP da);
void cpglm_fitted_x(double *x, SEXP da);

// cplm
double cplm_post_latT(double x, double y, double phi, double p) ; 
double cplm_lambda_tpois(double y, double phi, double p) ;
// function to simulate from 0 truncated poisson
int cplm_rtpois(double lambda ) ;
// function to compute the density of 0-truncated poisson on log scale
double cplm_dtpois(double x, double lambda) ;
void cplm_rlatT_reject(SEXP da) ;
double cplm_llik_lat(SEXP da);


// M-H 
int metrop_mvnorm_rw(int d, double *m, double *v, double *sn, 
		     double (*myfunc)(double *x, void *data), 
		     void *data);

int metrop_tnorm_rw( double m, double sd, double lb, double rb, double *sn, 
		     double (*myfunc)(double x, void *data), 
		     void *data) ;
void cov(int n, int p, double *x, double *ans) ;
void solve_po(int d, double *v, double *iv) ;
void mult_xtx(int m, int n, double *x, double *out) ;
double dmvnorm(int d, double *x, double *m, double *iv) ;

// tweedie 
double dtweedie(double y, double mu,
                 double phi, double p);
double dl2tweedie(int n, double *y, double *mu,
                  double phi, double p) ;
void rwishart(int d, double nu, double *scal, double *out) ;

// utility functions in bcpglmm
void bcpglmm_set_sims(SEXP da, int ns, double **sims) ;
void bcpglmm_set_init(SEXP da, int k) ;
void cpglmm_fitted(SEXP da) ;
void cpglmm_fitted_bx(double *x, SEXP da) ;
void cpglmm_fitted_ux(double *x, SEXP da) ;


// functions to export 
SEXP cpglmm_optimize(SEXP x) ;
SEXP cpglmm_update_mu(SEXP x) ;
SEXP cpglmm_update_u(SEXP x) ;
SEXP cpglmm_update_L(SEXP x) ;
SEXP cpglmm_update_dev(SEXP x, SEXP pm) ;
SEXP cpglmm_setPars(SEXP x, SEXP pm) ;
SEXP bcpglmm_gibbs_tw (SEXP da) ;
SEXP bcpglmm_gibbs_lat (SEXP da) ;
SEXP bcpglm_gibbs_lat (SEXP da) ;
SEXP bcpglm_gibbs_tw (SEXP da) ;

// numerical derivatives
void grad(int n, double *x, 
          double (*myfunc)(double *x, void *data), 
          void *data, double *ans) ;
void hess(int n, double *x, 
          double (*myfunc)(double *x, void *data), 
          void *data, double *ans) ;
