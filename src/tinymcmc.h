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
#include "stdlib.h"
#include "stdio.h"
#include "../randamu/include/rng.h"    

#define TM_MCMC_BASE 0  /*< The base Metropolis-Hastings algorithm */
#define TM_MCMC_PARALLEL_TEMPERING 1 /*< Like base, but with parallel tempering */

#define INVALID_NUMBER (NAN)    
#define IS_FINITE(x) (!(isnan(x) || isinf(x)))
#define TM_STOP -1 /* Stop the MCMC algorithm */    
#define TM_OK 0 /* No error */
#define TM_INFINITY_ENCOUNTERED 1 /* An infinity or NaN were encountered */
#define MIN(x, y) ((x) < (y) ? x : y)
#define MAX(x, y) ((x) > (y) ? x : y)    
    typedef int tm_error;
    typedef int tm_type;

    typedef struct {
        /**
         An array containing len elements, stored in v.
         */
        double* v;
        int len;
        int capacity;
    }
    tm_array;

    typedef struct {
        /**
         *  This contains the parameter values for each chain. 
         */
        int id;
        int length;
        int n;
        tm_array* values;
        tm_array** log_likelihood;
        tm_array* steps;
        tm_array* position;
        double merit;
        tm_array* next;
        double beta;
        rd_rng rng;
    } tm_chain;

    typedef struct {
    } tm_statistics;

    struct _tm_functions;
    typedef struct _tm_functions tm_functions;
    struct _tm_options;
    typedef struct _tm_options tm_options;

    typedef struct {
        /**
         * This struct contains the current state of MCMC.
         */
        int nchains;
        tm_chain** chains;
        tm_statistics* stats;
        tm_functions* functions;
        tm_options* options;
        void* user;
    } tm_state;

    struct _tm_prior;
    typedef struct _tm_prior tm_prior;

    typedef double (*tm_likelihood)(tm_array*, tm_state*);
    /*< The likelihood function. */
    typedef double (*tm_prior_function)(double, tm_prior*, tm_state*);
    /*< A prior function. */
    typedef void (*tm_proposal)(tm_array*, tm_array*, tm_state*, tm_chain*);
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
        tm_prior_function prior_function;
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
        tm_prior** priors;
        /**< Prior functions, specified for each parameter (required) */
        tm_proposal* proposal;
        /**< Proposal function */
        tm_callback* callback;
        /**< A callback. */
    };

    struct _tm_options {
        int chains;
        /**< Number of chains to use to check convergence on */
        double acceptance_rate_goal;
        /**< Goal acceptance rate (used to determine steps) */
        int burn_in;
        /**< Number of steps to discard at the beginning */
        int discard;
        /**< Number of steps to discard to minimize correlations */
        int exchanges;
        /**< For PT, how often to try an exchange */
    };


    /* Prior functions */
    tm_prior* tm_new_prior(tm_prior_function function, double min, double max, void* user);
    void tm_free_prior(tm_prior* prior);

    tm_prior* tm_uniform_prior(double min, double max);
    tm_prior* tm_modified_jeffreys_prior(double max, double a);
    tm_prior* tm_jeffreys_prior(double min, double max);

    /* Proposal functions */
    void tm_gaussian_gibbs1_proposal(tm_array* x0, tm_array* x1, tm_state* state, tm_chain* chain);

    /* Initialize algorithm */
    tm_state* tm_init(int type, tm_functions* functions, tm_options* options,
            tm_error* error);
    void tm_free(tm_state* state);

    /* Read/write state */
    tm_error tm_write(FILE* out, tm_state* state);
    tm_state* tm_read(FILE* in, tm_error* error);

    /* Actions */
    int tm_burn_in(tm_state* state, tm_error* error);
    int tm_run(tm_state* state, int steps, tm_error* error);

#ifdef	__cplusplus
}
#endif

#endif	/* TINYMCMC_H */

