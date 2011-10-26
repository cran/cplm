#include "cplm.h"
#include <R_ext/Rdynload.h>
#include "Matrix.h"

/** cholmod_common struct local to the cplm package */
cholmod_common c;

/** This is the CHOLMOD error handler from lme4*/
void attribute_hidden
cplm_R_cholmod_error(int status, const char *file, int line, const char *message)
{
    if(status < 0) {
#ifdef Failure_in_matching_Matrix
/* This fails unexpectedly with
 *  function 'cholmod_l_defaults' not provided by package 'Matrix'
 * from ../tests/lmer-1.R 's  (l.68)  lmer(y ~ habitat + (1|habitat*lagoon)
 */
	M_cholmod_defaults(&c);/* <--- restore defaults,
				* as we will not be able to .. */
	c.final_ll = 1;	    /* LL' form of simplicial factorization */
#endif

	error(_("Cholmod error '%s' at file:%s, line %d"), message, file, line);
    }
    else
	warning("Cholmod warning '%s' at file:%s, line %d",
		message, file, line);
}

/** Initializer for cplm, called upon loading the package.
 *
 *  Initialize CHOLMOD and require the LL' form of the factorization.
 *  Install the symbols to be used by functions in the package.
 */
#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
void R_init_cplm(DllInfo *dll)
{
    M_R_cholmod_start(&c);
    c.final_ll = 1;	    /* LL' form of simplicial factorization */

    /* need own error handler, that resets  final_ll (after *_defaults()) : */
    c.error_handler = cplm_R_cholmod_error;

}

/** Finalizer for lme4 called upon unloading the package.
 *
 */
void R_unload_cplm(DllInfo *dll){
    M_cholmod_finish(&c);
}
