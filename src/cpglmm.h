#ifndef CPGLMM_H
#define CPGLMM_H

/* cpglmm */
SEXP cpglmm_optimize(SEXP x) ;
SEXP cpglmm_update_mu(SEXP x) ;
SEXP cpglmm_update_u(SEXP x) ;
SEXP cpglmm_update_L(SEXP x) ;
SEXP cpglmm_update_dev(SEXP x, SEXP pm) ;
SEXP cpglmm_setPars(SEXP x, SEXP pm) ;

#endif
