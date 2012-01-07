/************************************************************/
/*   Function for fitting the Compound Poisson Generalized  */
/*    Linear Mixed models using Laplace approximation       */
/*    This is based the code from the R package lme4.       */
/*              Author:  Wayne Zhang                        */
/*            actuary_zhang@hotmail.com                     */
/************************************************************/


/**
 * @file cpglmm_lap.c
 * @brief Function for the Laplacian approximation in the        
 * Compound Poisson Generalized Linear Mixed Model.        
 * The program is based on the C code from lme4, where    
 * I removed the unnecessary pieces for linear mixed      
 * model and nonlinear model, and changed the code for    
 * updating the mean and the deviance and optimization     
 * @author Wayne Zhang                        
 * @date 2011-10-15    
*/

#include "cplm.h"
#include <R_ext/stats_package.h> /* for S_nlminb_iterate */
#include "Matrix.h"		 /* for cholmod functions */


/** cholmod_common struct initialized  */
extern cholmod_common c;

/* Constants used in Penalized Least Squares*/
/** Maximum number of iterations in update_u */
#define CM_MAXITER  300
/** Tolerance level for convergence criterion in update_u */
#define CM_TOL      1e-10
/** Minimum step factor in update_u */
#define CM_SMIN     1e-5


/** positions in the deviance vector */
enum devP {
    ML_POS=0,			/**<Maximum likelihood estimation criterion  */
    REML_POS,			/**<REML criterion */
    ldL2_POS,			/**<2*log-determinant of L */
    ldRX2_POS,			/**<2*log-determinant of RX */
    sigmaML_POS,		/**<current ML estimate of sigma */
    sigmaREML_POS,		/**<current REML estimate of sigma */
    pwrss_POS,			/**<penalized weighted residual sum of squares */
    disc_POS,			/**<discrepancy */
    usqr_POS,			/**<squared length of u */
    wrss_POS,			/**<weighted residual sum of squares  */
    dev_POS,			/**<deviance - defined for quasi families  */
    llik_POS,			/**<log-likelihood - undefined for quasi families  */
    NULLdev_POS			/**<null deviance */
};

/** positions in the dims vector */
enum dimP {
    nt_POS=0,			/**<number of terms in random effects */
    n_POS,			/**<number of observations */
    p_POS,			/**<number of fixed-effects parameters */
    q_POS,			/**<number of random effects */
    s_POS,			/**<number of variables in h (1 unless nonlinear) */
    np_POS,			/**<total number of parameters for T and S */
    LMM_POS,			/**<is the model a linear mixed model? */
    isREML_POS,			/**<indicator of REML estimation */
    fTyp_POS,			/**<family type for generalized model */
    lTyp_POS,			/**<link type for generalized model */
    vTyp_POS,			/**<variance type for generalized model */
    nest_POS,			/**<indicator of nested grouping factors */
    useSc_POS,			/**<does the family use a separate scale parameter */
    nAGQ_POS,			/**<number of adaptive Gauss-Hermite quadrature pts */
    verb_POS,			/**<verbose output in mer_optimize? */
    mxit_POS,			/**<maximum # of iterations in mer_optimize */
    mxfn_POS,			/**<maximum # of function evaluations in mer_optimize */
    cvg_POS			/**<convergence indictor from port optimization  */
};

/**
 * Extract the slot named nm from the object obj and return a null pointer
 * if the slot has length 0 or a pointer to the REAL contents.
 *
 * @param obj pointer to an S4 object
 * @param str pointer to a symbol naming the slot to extract
 *
 * @return pointer to the REAL contents, if nonzero length, otherwise
 * a NULL pointer
 *
 */
static R_INLINE double *SLOT_REAL_NULL(SEXP obj, char *str)
{
    SEXP pt = GET_SLOT(obj, install(str));
    return LENGTH(pt) ? REAL(pt) : (double*) NULL; 
}

/** Return the double pointer to the X slot */
#define X_SLOT(x) SLOT_REAL_NULL(x, "X")

/** Return the double pointer to the y slot */
#define Y_SLOT(x) SLOT_REAL_NULL(x, "y")

/** Allocate (alloca) a cholmod_sparse struct, populate it with values
 * from the Zt slot and return the pointer. */
#define Zt_SLOT(x) AS_CHM_SP(GET_SLOT(x, install("Zt")))

/** Return the double pointer to the offset slot or (double*) NULL if
 * offset has length 0) */
#define OFFSET_SLOT(x) SLOT_REAL_NULL(x, "offset")

/** Return the double pointer to the pWt slot or (double*) NULL if
 * pWt has length 0) */
#define PWT_SLOT(x) SLOT_REAL_NULL(x, "pWt")

/** Return the integer pointer to the dims slot */
#define DIMS_SLOT(x) INTEGER(GET_SLOT(x, install("dims")))

/** Return the double pointer to the fixef slot */
#define FIXEF_SLOT(x) SLOT_REAL_NULL(x, "fixef")

/** Return the double pointer to the u slot */
#define U_SLOT(x) SLOT_REAL_NULL(x, "u")

/** Return the double pointer to the eta slot */
#define ETA_SLOT(x) SLOT_REAL_NULL(x, "eta")

/** Return the double pointer to the mu slot */
#define MU_SLOT(x) SLOT_REAL_NULL(x, "mu")

/** Return the double pointer to the muEta slot or (double*) NULL if
 * muEta has length 0) */
#define MUETA_SLOT(x) SLOT_REAL_NULL(x, "muEta")

/** Return the double pointer to the var slot or (double*) NULL if
 * var has length 0) */
#define VAR_SLOT(x) SLOT_REAL_NULL(x, "var")

/** Return the double pointer to the resid slot */
#define RESID_SLOT(x) SLOT_REAL_NULL(x, "resid")

/** Allocate (alloca) a cholmod_sparse struct, populate it with values
 * from the A slot and return the pointer. */
#define A_SLOT(x) AS_CHM_SP(GET_SLOT(x, install("A")))

/** Allocate (alloca) a cholmod_factor struct, populate it with values
 * from the L slot and return the pointer. */
#define L_SLOT(x) AS_CHM_FR(GET_SLOT(x, install("L")))

/** Return the integer pointer to the Gp slot */
#define Gp_SLOT(x) INTEGER(GET_SLOT(x, install("Gp")))

/** Return the double pointer to the Cx slot or (double*) NULL if
 * Cx has length 0) */
#define Cx_SLOT(x) SLOT_REAL_NULL(x, "Cx")

/** Return the double pointer to the deviance slot */
#define DEV_SLOT(x) SLOT_REAL_NULL(x, "deviance")

/** Return the double pointer to the sqrtrWt slot or (double*) NULL if
 *  sqrtrWt has length 0) */
#define SRWT_SLOT(x) SLOT_REAL_NULL(x, "sqrtrWt")

/** Return the double pointer to the sqrtXWt slot or (double*) NULL if
 *  sqrtXWt has length 0) */
#define SXWT_SLOT(x) SLOT_REAL_NULL(x, "sqrtXWt")

/** Return the double pointer to the p slot or (double*) NULL if
 *  sqrtXWt has length 0) */
#define P_SLOT(x) SLOT_REAL_NULL(x, "p")

/** Return the double pointer to the phi slot or (double*) NULL if
 *  sqrtXWt has length 0) */
#define PHI_SLOT(x) SLOT_REAL_NULL(x, "phi")

/** Return the double pointer to the link_power slot or (double*) NULL if
 *  sqrtXWt has length 0) */
#define LKP_SLOT(x) SLOT_REAL_NULL(x, "link.power")

/** Return the double pointer to the bound_p slot  */
#define BDP_SLOT(x) SLOT_REAL_NULL(x,"bound.p")

/** Return the integer pointer to the permutation vector in the L slot */
#define PERM_VEC(x) INTEGER(GET_SLOT(GET_SLOT(x, install("L")), install("perm")))

/** Return the double pointer to the ranef slot or (double*) NULL if
 *  ranef has length 0) */
#define RANEF_SLOT(x) SLOT_REAL_NULL(x, "ranef")

/** Return the double pointer to the ghw slot */
#define GHW_SLOT(x) SLOT_REAL_NULL(x, "ghw")

/** Return the double pointer to the ghx slot */
#define GHX_SLOT(x) SLOT_REAL_NULL(x, "ghx")

/**
 * Permute the vector src according to perm into dest
 *
 * @param dest destination
 * @param src source
 * @param perm NULL or 0-based permutation of length n
 * @param n length of src, dest and perm
 *
 * @return dest
 *
 * \note If perm is NULL the first n elements of src are copied to dest.
 */
static R_INLINE double*
apply_perm(double *dest, const double *src, const int *perm, int n)
{
    for (int i = 0; i < n; i++) dest[i] = src[perm ? perm[i] : i];
    return dest;
}

/**
 * Return the sum of squares of the first n elements of x
 *
 * @param n
 * @param x
 *
 * @return sum of squares
 */
static R_INLINE double sqr_length(const double *x, int n)
{
    double ans = 0;
    for (int i = 0; i < n; i++) ans += x[i] * x[i];
    return ans;
}

/**
 * Return the index of the term associated with parameter index ind
 *
 * @param ind an index in [0, Gp[nt] - 1]
 * @param nt total number of terms
 * @param Gp group pointers, a vector of length nt+1 with Gp[0] = 0
 *
 * @return sum of squares
 */
static R_INLINE int Gp_grp(int ind, int nt, const int *Gp)
{
    for (int i = 0; i < nt; i++) if (ind < Gp[i + 1]) return i;
    error(_("invalid row index %d (max is %d)"), ind, Gp[nt]);
    return -1;                  /* -Wall */
}


/**
 * Extract the parameters from ST list
 *
 * @param x a list object
 * @param pars vector of the appropriate length
 *
 * @return pointer to the parameter vector
 */
static double *ST_getPars(SEXP x, double *pars)
{
    SEXP ST = GET_SLOT(x, install("ST"));
    int nt = LENGTH(ST), pos = 0;
    for (int i = 0; i < nt; i++) {
	SEXP STi = VECTOR_ELT(ST, i);
	double *st = REAL(STi);
	int nci = INTEGER(getAttrib(STi, R_DimSymbol))[0];
	int ncp1 = nci + 1;

	for (int j = 0; j < nci; j++)
	    pars[pos++] = st[j * ncp1];  /**<get parm along the diagonal (T) */
	for (int j = 0; j < (nci - 1); j++) //get S
	    for (int k = j + 1; k < nci; k++)
		pars[pos++] = st[k + j * nci];
    }
    return pars;
}

/**
 * Populate the st, nc and nlev arrays.  Return the maximum element of nc.
 *
 * @param ST pointer to a list (length nt) of matrices
 * @param Gp group pointers (length nt + 1)
 * @param st length nt array of (double*) pointers to be filled with
 * pointers to the contents of the matrices in ST.  Not used if NULL.
 * @param nc length nt array to be filled with the number of columns
 * @param nlev length nt array to be filled with the number of
 *        levels of the grouping factor for each term
 *
 * @return maximum element of nc
 */
static int			/* populate the st, nc and nlev arrays */
ST_nc_nlev(const SEXP ST, const int *Gp, double **st, int *nc, int *nlev)
{
    int ans = 0, nt = LENGTH(ST);

    for (int i = 0; i < nt; i++) {
	SEXP STi = VECTOR_ELT(ST, i);
	int nci = *INTEGER(getAttrib(STi, R_DimSymbol));

	if (nci > ans) ans = nci;
	if (st) st[i] = REAL(STi);
	nc[i] = nci;
	nlev[i] = (Gp[i + 1] - Gp[i])/nci;
    }
    return ans;
}

/**
 * Multiply A on the left by T'
 *
 * @param A sparse model matrix
 * @param Gp group pointers
 * @param nc number of columns per term
 * @param nlev number of levels per term
 * @param st ST arrays for each term
 * @param nt number of terms
 *
 */
static void Tt_Zt(CHM_SP A, int *Gp, int *nc, int *nlev, double **st, int nt)
{
    int *ai = (int*)(A->i), *ap = (int *)(A->p);
    double *ax = (double*)(A->x), one[] = {1,0};

    for (int j = 0; j < A->ncol; j++) /* multiply column j by T' */
	for (int p = ap[j]; p < ap[j + 1];) {
	    int i = Gp_grp(ai[p], nt, Gp);

	    if (nc[i] <= 1) p++;
	    else {
		int nr = p;	/* number of rows in `B' in dtrmm call */
		while ((ai[nr] - Gp[i]) < nlev[i]) nr++;
		nr -= p;	/* nr == 1 except in models with carry-over */
		F77_CALL(dtrmm)("R", "L", "N", "U", &nr, nc + i,
				one, st[i], nc + i, ax + p, &nr);
		p += (nr * nc[i]);
	    }
	}
}

/**
 * Evaluate the sparse model matrix A from the Zt and ST slots
 * A = S'T'Z', the transpose of the design matrix with unit random effects 
 *
 * @param x a list object
 */
static void update_A(SEXP x)
{
    CHM_SP A = A_SLOT(x), Zt = Zt_SLOT(x);
    int *Gp = Gp_SLOT(x), *ai = (int*)(A->i), *ap = (int*)(A->p),
	 *zp = (int*)(Zt->p), ncmax,
	nt = DIMS_SLOT(x)[nt_POS];
    int annz = ap[A->ncol], znnz = zp[Zt->ncol];
    int *nc = Alloca(nt, int), *nlev = Alloca(nt, int);
    double **st = Alloca(nt, double*), *ax = (double*)(A->x),
	*zx = (double*)(Zt->x);
    R_CheckStack();

    ncmax = ST_nc_nlev(GET_SLOT(x, install("ST")), Gp, st, nc, nlev);

     /* Copy Z' to A unless A has new nonzeros */
    Memcpy(ax, zx, znnz);
    /* When T != I multiply A on the left by T' */
    if (ncmax > 1) Tt_Zt(A, Gp, nc, nlev, st, nt);
				/* Multiply A on the left by S */
    for (int p = 0; p < annz; p++) {
	int i = Gp_grp(ai[p], nt, Gp);
	ax[p] *= st[i][((ai[p] - Gp[i]) / nlev[i]) * (nc[i] + 1)];
    }
}


/**
 * Update the L, sqrtrWt and resid slots.  It is assumed that
 * update_mu has already been called at the current values of u and
 * the model parameters and that A has been updated.
 *
 * @param x pointer to an mer object
 *
 * @return penalized weighted residual sum of squares
 *
 */
static double cp_update_L(SEXP x)
{
    int *dims = DIMS_SLOT(x);
    int n = dims[n_POS] ;
    double  *cx = Cx_SLOT(x), *d = DEV_SLOT(x),
	*res = RESID_SLOT(x), *mu = MU_SLOT(x), *muEta = MUETA_SLOT(x),
	*pwt = PWT_SLOT(x), *sXwt = SXWT_SLOT(x), *srwt = SRWT_SLOT(x),
	*var =  VAR_SLOT(x), *y = Y_SLOT(x), one[] = {1,0};
    CHM_SP A = A_SLOT(x);
    CHM_FR L = L_SLOT(x);
    R_CheckStack();

    // Update srwt and res. Reevaluate wrss. 
    d[wrss_POS] = 0;
    for (int j = 0; j < n; j++) {
        srwt[j] = sqrt((pwt ? pwt[j] : 1.0) / (var ? var[j] : 1.0)) ;
        res[j] = srwt[j] * (y[j] - mu[j]);
        d[wrss_POS] += res[j] * res[j];
    }
    int  *ap = (int*)A->p;
    double *ax = (double*)(A->x);

    // sXwt is the sqrt of the weight in the penalized regression
    for (int i = 0; i < n; i++) 
          sXwt[i] = (srwt ? srwt[i] : 1) * (muEta ? muEta[i] : 1) ;
    
    //  a scaled version of A: A * sqrt(W) 
    for (int j = 0; j < n; j++)
        for (int p = ap[j]; p < ap[j + 1]; p++)
            cx[p] = ax[p] * sXwt[j] ;
    A->x = (void*)cx;
    
    // compute the cholesky factor of AA'+I with permutation
    if (!M_cholmod_factorize_p(A, one, (int*)NULL, 0, L, &c))
	error(_("cholmod_factorize_p failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);
    
    d[ldL2_POS] = M_chm_factor_ldetL2(L);
    d[pwrss_POS] = d[usqr_POS] + d[wrss_POS];
    return d[pwrss_POS];
}


/**
 * Update the eta, v, mu, resid and var slots according to the current
 * values of the parameters and u.  Also evaluate d[wrss_POS] using
 * the current sqrtrWt slot.  The sqrtrWt slot is changed in update_L.
 * This function is different from that in lme4, and hence the prefix
 * "cp" is added. Specifically, the update of mu, muEta and variance
 * is from a Tweedie compound Poisson model.
 *
 * @param x pointer to a list object
 *
 * @return penalized, weighted residual sum of squares
 */
static double cp_update_mu(SEXP x)
{
    int *dims = DIMS_SLOT(x);
    int i1 = 1, n = dims[n_POS], p = dims[p_POS];
    int  i ;
    double  *d = DEV_SLOT(x), *eta = ETA_SLOT(x),
        *mu = MU_SLOT(x),
	*muEta = MUETA_SLOT(x), *offset = OFFSET_SLOT(x),
	*srwt = SRWT_SLOT(x), *res = RESID_SLOT(x),
        *lkp = LKP_SLOT(x), *vp = P_SLOT(x),
        *var = VAR_SLOT(x), 
        *y = Y_SLOT(x), one[] = {1,0};
    CHM_FR L = L_SLOT(x);
    CHM_SP A = A_SLOT(x);
    CHM_DN Ptu, ceta, cu = AS_CHM_DN(GET_SLOT(x, install("u")));
    R_CheckStack();

    /* eta := offset or eta := 0 */
    for ( i = 0; i < n; i++) eta[i] = offset ? offset[i] : 0;
				/* eta := eta + X \beta */
    F77_CALL(dgemv)("N", &n, &p, one, X_SLOT(x), &n,
		    FIXEF_SLOT(x), &i1, one, eta, &i1);
				/* eta := eta + C' P' u */
    Ptu = M_cholmod_solve(CHOLMOD_Pt, L, cu, &c);
    ceta = N_AS_CHM_DN(eta, n, 1);
    R_CheckStack();
    if (!M_cholmod_sdmult(A, 1 , one, one, Ptu, ceta, &c))
	error(_("cholmod_sdmult error returned"));
    M_cholmod_free_dense(&Ptu, &c);

    // update mu, muEta and var
    cplm_mu_eta(mu, muEta, n, eta, *lkp);
    cplm_varFun(var,mu, n,*vp);
   
    d[wrss_POS] = 0;		/* update resid slot and d[wrss_POS] */
    for (i = 0; i < n; i++) {
	res[i] = (y[i] - mu[i]) * (srwt ? srwt[i] : 1);
	d[wrss_POS] += res[i] * res[i];
    }
				/* store u'u */
    d[usqr_POS] = sqr_length((double*)(cu->x), dims[q_POS]);
    d[pwrss_POS] = d[usqr_POS] + d[wrss_POS];
    return d[pwrss_POS];
}



/**
 * Iterate to determine the conditional modes of the random effects.
 * This is the same as that in lme4. 
 * @param x pointer to a list object
 *
 * @return number of iterations to convergence (0 for non-convergence)
 */

static int cp_update_u(SEXP x)
{
    int *dims = DIMS_SLOT(x);
    int i, n = dims[n_POS], q = dims[q_POS], verb = dims[verb_POS];
    double *Cx = Cx_SLOT(x), 
	*res = RESID_SLOT(x), *u = U_SLOT(x),
	cfac = ((double)n) / ((double)q),
	crit, pwrss, pwrss_old, step;
    double *tmp = Calloc(q, double), *tmp1 = Calloc(q, double),
	*uold = Calloc(q, double), one[] = {1,0}, zero[] = {0,0};
    CHM_FR L = L_SLOT(x);
    CHM_DN cres = N_AS_CHM_DN(res, n, 1),
	ctmp = N_AS_CHM_DN(tmp, q, 1), sol;
    
    R_CheckStack();

    if (!(L->is_ll)) error(_("L must be LL', not LDL'"));
    CHM_SP C  = A_SLOT(x);
    R_CheckStack();
    C->x = (void*)Cx;
    
    // update mu related parms
    AZERO(u, q);
    cp_update_mu(x);
 
    for (i = 0; ; i++) {
	Memcpy(uold, u, q);
	pwrss_old = cp_update_L(x);
        // tmp := PC %*% wtdResid 
	M_cholmod_sdmult(C, 0 , one, zero, cres, ctmp, &c);
	Memcpy(tmp1, tmp, q);
	apply_perm(tmp, tmp1, (int*)L->Perm, q);
				// tmp := tmp - u 
	for (int j = 0; j < q; j++) tmp[j] -= u[j];
				// solve L %*% sol = tmp 
	if (!(sol = M_cholmod_solve(CHOLMOD_L, L, ctmp, &c)))
	    error(_("cholmod_solve (CHOLMOD_L) failed"));
	Memcpy(tmp, (double*)(sol->x), q);
	M_cholmod_free_dense(&sol, &c);
				// evaluate convergence criterion 
	crit = cfac * sqr_length(tmp, q) / pwrss_old;
	if (crit < CM_TOL) break; // don't do needless evaluations 
				// solve t(L) %*% sol = tmp 
	if (!(sol = M_cholmod_solve(CHOLMOD_Lt, L, ctmp, &c)))
	    error(_("cholmod_solve (CHOLMOD_Lt) failed"));
	Memcpy(tmp, (double*)(sol->x), q);
	M_cholmod_free_dense(&sol, &c);
	for (step = 1; step > CM_SMIN; step /= 2) { // step halving 
	    for (int j = 0; j < q; j++) u[j] = uold[j] + step * tmp[j];
	    pwrss = cp_update_mu(x);
	    if (verb < 0)
		Rprintf("%2d,%8.6f,%12.4g: %15.6g %15.6g %15.6g %15.6g\n",
			i, step, crit, pwrss, pwrss_old, u[1], u[2]);
	    if (pwrss < pwrss_old) {
		pwrss_old = pwrss;
		break;
	    }
	}
	if (step <= CM_SMIN || i > CM_MAXITER) return 0;
    }
    
    Free(tmp); Free(tmp1); Free(uold);
    return i;
}

/**
 * Update the ST and A slots of an mer object.
 *
 * @param x an mer object
 * @param pars double vector of the appropriate length
 * @return updated deviance
 *
 */
static void
ST_setPars(SEXP x, const double *pars)
{
    int *Gp = Gp_SLOT(x), nt = DIMS_SLOT(x)[nt_POS], pos = 0;
    int *nc = Alloca(nt, int), *nlev = Alloca(nt, int);
    double **st = Alloca(nt, double*);
    R_CheckStack();

    ST_nc_nlev(GET_SLOT(x, install("ST")), Gp, st, nc, nlev);
				/* install the parameters in the ST slot */
    for (int i = 0; i < nt; i++) {
	int nci = nc[i], ncp1 = nc[i] + 1;
	double *sti = st[i];

	for (int j = 0; j < nci; j++)
	    sti[j * ncp1] = pars[pos++];
	for (int j = 0; j < (nci - 1); j++)
	    for (int k = j + 1; k < nci; k++)
		sti[k + j * nci] = pars[pos++];
    }
    update_A(x);
}


/**
 * Update the ST, A, fixef, phi and p slots 
 *
 * @param x an cpglmm object
 * @param pars double vector of the appropriate length: theta, beta, log(phi) and p in order
 * @return updated deviance
 *
 */
static void cp_setPars(SEXP x, double *pars){
    int *dims = DIMS_SLOT(x) ;
    int phi_POS = dims[np_POS] + dims[p_POS] ;    
    double  *phi = PHI_SLOT(x),
      *p = P_SLOT(x), *fixef = FIXEF_SLOT(x) ;
    // update ST from parm
    ST_setPars(x, pars);
    // update fixef, phi and p
    Memcpy(fixef, pars + dims[np_POS], dims[p_POS]);
    phi[0] = exp(pars[phi_POS])  ;
    p[0] = pars[phi_POS + 1] ;       
}

/**
 * Compute twice negative marginal logliklihood to be used in optimization. 
 * This is different from that in lme4, so I added the prefix "cp".
 *
 * @param x a cpglmm object
 * @param parm vector of parameters: theta, beta, log(phi) and p in order
 *            
 * @return twice negative marginal logliklihood
 *
 */
static double cp_update_dev(SEXP x, double *parm)
{
    int *dims = DIMS_SLOT(x) ;
    int n = dims[n_POS], 
      q = dims[q_POS], nAGQ = dims[nAGQ_POS] ;
    SEXP flistP = GET_SLOT(x, install("flist"));
    double *d = DEV_SLOT(x), *y = Y_SLOT(x),
      *mu = MU_SLOT(x), *phi = PHI_SLOT(x),
      *p = P_SLOT(x), *pwt = PWT_SLOT(x),
      *u = U_SLOT(x);
    CHM_FR L = L_SLOT(x);
    
    if (parm != NULL) cp_setPars(x, parm) ;
    // find conditional mode
    cp_update_u(x);

    if (nAGQ < 1) error("nAGQ must be positive");
    if ((nAGQ > 1) & (LENGTH(flistP) != 1))
      error("AGQ method requires a single grouping factor");
    d[ML_POS] = d[ldL2_POS];

    // Laplace Approximation 
    if (nAGQ == 1) {
      double ans=0 ;
      for (int i=0;i<n;i++)
        ans += dtweedie(y[i],mu[i],phi[0]/pwt[i],p[0]) ; //log-likelihood
      d[disc_POS] = - 2* ans;
      d[ML_POS] += d[disc_POS] + d[usqr_POS] / phi[0] ; 
      return d[ML_POS]; 
    } else {
      // Adaptive Gauss-Hermite quadrature 
      const int nl = nlevels(VECTOR_ELT(flistP, 0));
      const int nre = q / nl; 	/* number of random effects per level */
      int *fl0 = INTEGER(VECTOR_ELT(flistP, 0)), *pointer = Alloca(nre, int);
      double *ghw = GHW_SLOT(x), *ghx = GHX_SLOT(x),
	//
	*uold = Memcpy(Calloc(q, double), u, q),
	*tmp = Calloc(nl, double),
	w_pro = 1, z_sum = 0;           /* values needed in AGQ evaluation */
      const double sigma = sqrt(phi[0]);
      R_CheckStack();
      AZERO(pointer, nre);
      AZERO(tmp, nl);

      while(pointer[nre - 1] < nAGQ){
	double *z = Calloc(q, double);       /* current abscissas */
	double *ans = Calloc(nl, double);    /* current penalized residuals in different levels */
	/* update abscissas and weights */
	for(int i = 0; i < nre; ++i){
	  for(int j = 0; j < nl; ++j)
	    z[i + j * nre] = ghx[pointer[i]];
	  w_pro *= ghw[pointer[i]];
	  z_sum += z[pointer[i]] * z[pointer[i]]; // constant used in AGQ
	}
	// scaling: u = u_old + sigma * z / L_ii
	CHM_DN cz = N_AS_CHM_DN(z, q, 1), sol;
	if(!(sol = M_cholmod_solve(CHOLMOD_L, L, cz, &c)))
	  error(_("cholmod_solve(CHOLMOD_L) failed"));
	Memcpy(z, (double *)sol->x, q);
	M_cholmod_free_dense(&sol, &c);   // z / L_ii
	for(int i = 0; i < q; ++i) 
	  u[i] = uold[i] + sqrt(2) * sigma * z[i];   
	cp_update_mu(x);
	// get data likelihood
	AZERO(ans, nl);
	for (int i = 0; i < n; i++) {
	  ans[fl0[i] - 1] +=  dtweedie(y[i], mu[i], phi[0] / pwt[i], p[0]) ; 
	}
	// contribution from u
	for(int i = 0; i < nre; ++i)
	  for(int j = 0; j < nl; ++j)
	    ans[j] -= 0.5 * u[i + j * nre] * u[i + j * nre] / phi[0];
	// adaptive GHQ estimate 
	for(int i = 0; i < nl; ++i)
	  tmp[i] += exp(ans[i] + z_sum ) * w_pro / sqrt(PI);
	// move pointer to next combination of weights and abbsicas 
	int count = 0;
	pointer[count]++;
	while(pointer[count] == nAGQ && count < nre - 1){
	  pointer[count] = 0;
	  pointer[++count]++;
	}
	w_pro = 1;
	z_sum = 0;
	if(z)    Free(z);
	if(ans)  Free(ans);
      }
      for(int j = 0; j < nl; ++j) d[ML_POS] -= 2 * log(tmp[j]);
      Memcpy(u, uold, q);
      cp_update_mu(x);
      if(tmp)   Free(tmp);
      if(uold)  Free(uold);
      return d[ML_POS] ;
    }
}


/**
 * Optimize the  Laplace approximation to the deviance
 * of a cpglmm object. The parameter set includes
 * theta (variance component), beta (fixed effects),
 * phi (scale parameter) and p (index parameter) 
 *
 * @param x pointer to an mer object
 *
 * @return R_NilValue
 */
SEXP cpglmm_optimize(SEXP x)
{
    SEXP ST = GET_SLOT(x, install("ST"));
    int *dims = DIMS_SLOT(x);
    int  nt = dims[nt_POS];
    int nv1 = dims[np_POS] +  dims[p_POS], verb = dims[verb_POS];
    int nv = nv1 +2 ;
    int liv = S_iv_length(OPT, nv), lv = S_v_length(OPT, nv);
    double *g = (double*)NULL, *h = (double*)NULL, fx = R_PosInf;
    double *fixef = FIXEF_SLOT(x);
    int *iv = Alloca(liv, int);
    double *b = Alloca(2 * nv, double), *d = Alloca(nv, double),
	*v = Alloca(lv, double), *xv = Alloca(nv, double);
    R_CheckStack();

    // retrieve parameter values from x 
    ST_getPars(x, xv); /* update xv for theta */
    Memcpy(xv + dims[np_POS], fixef, dims[p_POS]); /*update xv for beta*/
    xv[nv1] = log(PHI_SLOT(x)[0]) ;
    xv[nv1+1] = P_SLOT(x)[0] ;
    
    double eta = 1.e-5; /* estimated rel. error on computed lpdisc */
    v[31] = eta;		/* RFCTOL */
    v[36] = eta;		/* SCTOL */
    v[41] = eta;		/* ETA0 */
				/* initialize the state vectors v and iv */
    S_Rf_divset(OPT, iv, liv, lv, v);
    iv[OUTLEV] = (verb < 0) ? -verb : verb;
    iv[MXFCAL] = dims[mxfn_POS];
    iv[MXITER] = dims[mxit_POS];
				/* set the bounds to plus/minus Infty  */
    for (int i = 0; i < nv; i++) {
	b[2*i] = R_NegInf; b[2*i+1] = R_PosInf; d[i] = 1;
    }
				/* reset lower bounds on elements of theta_S */
    for (int i = 0, pos = 0; i < nt; i++) {
	int nc = *INTEGER(getAttrib(VECTOR_ELT(ST, i), R_DimSymbol));
	for (int j = 0; j < nc; j++) b[pos + 2*j] = 0;
	pos += nc * (nc + 1);
    }
    /* reset bound for p */
    b[(nv1+1)*2]= BDP_SLOT(x)[0] ;
    b[nv*2-1] = BDP_SLOT(x)[1] ;

    /* run optimization */
    do {
        fx = cp_update_dev(x, xv);        
	S_nlminb_iterate(b, d, fx, g, h, iv, liv, lv, nv, v, xv);
    } while (iv[0] == 1 || iv[0] == 2);
    fx = cp_update_dev(x, xv) ; //update slots using the values at the minima
    dims[cvg_POS] = iv[0];
    return R_NilValue;
}


/* functions callable from R */

/**
 * R callable function to update L
 *
 * @param x a R list object
 *
 * @return penalized, weighted residual sum of squares
*/
SEXP cpglmm_update_L(SEXP x){
    return ScalarReal(cp_update_L(x));
}

/**
 * R callable function to update mu
 *
 * @param x a R list object
 *
 * @return penalized, weighted residual sum of squares
*/
SEXP cpglmm_update_mu(SEXP x){
    return ScalarReal(cp_update_mu(x));
}


/**
 * R callable function to update u
 *
 * @param x a R list object
 *
 * @return number of iterations to convergence (0 for non-convergence)
*/
SEXP cpglmm_update_u(SEXP x){
    return ScalarInteger(cp_update_u(x));
}

/**
 * R callable function to update deviance
 *
 * @param x a cpglmm object
 * @param pm vector of parameters: theta, beta, log(phi) and p in order
 *
 * @return twice negative loglikelihood
*/

SEXP cpglmm_update_dev(SEXP x, SEXP pm){  
  return ScalarReal(cp_update_dev(x, (pm != R_NilValue) ? REAL(pm) : (double*) NULL));            
}


/**
 * R callable function to update parameters
 *
 * @param x a cpglmm object
 * @param pm vector of parameters: theta, beta, log(phi) and p in order
 *
*/
SEXP cpglmm_setPars(SEXP x, SEXP pm){  
    cp_setPars(x, REAL(pm)) ;
    return R_NilValue ;      
}
