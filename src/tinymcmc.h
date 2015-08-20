/* 
 * File:   tinymcmc.h
 * Author: Stefano Meschiari
 * 
 * An inefficient, yet super simple implementation of the Metropolis-Hasting
 * algorithm.
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

    typedef struct {
        /**
         An array containing len elements, stored in v.
         */
        double* v;
        int len;
    } tm_array;

    typedef struct {
        /**
         *  This contains the parameter values for each chain. 
         */
        int id;
        int length;
        int capacity;
        double** values;
        double* log_likelihood;
        tm_array* steps;
        tm_array* position;
    } tm_chain;

    typedef struct {
    } tm_statistics;

    struct _tm_functions;
    typedef struct _tm_functions tm_functions;

    typedef struct {
        /**
         * This struct contains the current state of MCMC.
         */
        int nchains;
        tm_chain* chains;
        tm_statistics* stats;
        tm_functions* functions;
        void* user;
    } tm_state;

    struct _tm_prior;
    typedef struct _tm_prior tm_prior;

    typedef double (*tm_likelihood)(tm_array*, tm_chain*, void*);
    /*< The likelihood function. */
    typedef double (*tm_function)(double, tm_prior*, tm_state*);
    /*< A generic function. */
    typedef double (*tm_proposal)(tm_array*, tm_array*, tm_state*);
    /*< A proposal function. */
    typedef int (*tm_callback)(tm_state*);

    /*< A callback function. */

    struct _tm_prior {
        /**
         * A prior function, with limits [min:max] (inclusive)
         */
        double min;
        double max;
        void* user;
        tm_function prior_function;
        double* extra;
    };

    struct _tm_functions {
        /**
         * This struct contains the core functions needed by
         * the MCMC algorithm.
         */
        tm_array* parameters;
        /**< Starting values for parameters; pass NULL to randomly initialize  */
        tm_likelihood log_likelihood;
        /**< The likelihood function (required) */
        tm_array* steps;
        /**< Initial steps for each parameter; pass NULL to automatically calculate,
         * based on goal acceptance rate */
        tm_prior* priors;
        /**< Prior functions, specified for each parameter (required) */
    };

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


    tm_prior* tm_new_prior(tm_function function, double min, double max, void* user);
    void tm_free_prior(tm_prior* prior);

    tm_prior* tm_uniform_prior(double min, double max);
    tm_prior* tm_modified_jeffreys_prior(double max, double a);
    tm_prior* tm_jeffreys_prior(double min, double max);
#ifdef	__cplusplus
}
#endif

#endif	/* TINYMCMC_H */

