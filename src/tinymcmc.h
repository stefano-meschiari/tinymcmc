/* 
 * File:   tinymcmc.h
 * Author: sm52286
 *
 * Created on August 18, 2015, 1:34 PM
 */

#ifndef TINYMCMC_H
#define	TINYMCMC_H

#ifdef	__cplusplus
extern "C" {
#endif

#define TM_MCMC_BASE 0  /*< The base Metropolis-Hastings algorithm */
#define TM_MCMC_PARALLEL_TEMPERING 1 /*< Like base, but with parallel tempering */



    typedef double (*tm_likelihood)(double*, void*);
    typedef double (*tm_function)(double, void*);

    typedef struct {
        double* v;
        int len;
    } tm_array;

    typedef struct {
        double min;
        double max;
        tm_function prior_function;
    } tm_prior;

    typedef struct {
        /**
         * This struct contains the core functions needed by
         * the MCMC algorithm.
         */
        tm_array parameters;
        /**< Initial steps for each parameter */
        tm_likelihood log_likelihood;
        /**< The likelihood function (required) */
        tm_array steps;
        /**< Starting values for parameters; pass NULL to  */
        tm_prior* log_priors;
        /**< Prior functions, specified for each parameter (required) */
        void* user;
    }
    tm_functions;

    typedef struct {
        int chains;
        /**< Number of chains to use to check convergence on */
        double acceptance_rate_goal;
        /**< Goal acceptance rate (used to determine steps) */
        int burn_in;
        /**< Number of steps to discard at the beginning */
        int discard;
        /**< Number of steps to discard to minimize correlations */
    } tm_base_options;




#ifdef	__cplusplus
}
#endif

#endif	/* TINYMCMC_H */

