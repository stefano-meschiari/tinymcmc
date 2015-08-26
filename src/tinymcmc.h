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

#define TM_VERSION 0.001    

#define INVALID_NUMBER (NAN)    
#define IS_FINITE(x) (!(isnan(x) || isinf(x)))
#define TM_STOP -1 /* Stop the MCMC algorithm */    
#define TM_OK 0 /* No error */
#define TM_INFINITY_ENCOUNTERED 1 /* An infinity or NaN were encountered */
#define TM_COULD_NOT_WRITE_FILE 2 /* Could not open file for output */
#define TM_CRITICAL_ERROR 20
#define MIN(x, y) ((x) < (y) ? x : y)
#define MAX(x, y) ((x) > (y) ? x : y)    

#define TM_STAGE_RUN 0
#define TM_STAGE_BURNIN 1
#define TM_STAGE_STEPS 2    

    typedef int tm_error;
    typedef int tm_type;

    typedef struct {
        char str[8192];
    } tm_string;

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
        tm_array* x;
        double log_likelihood;
        double log_posterior;
    } tm_position;

    typedef struct {
        /**
         *  This contains the parameter values for each chain. 
         */
        int id;
        int length;
        tm_array** values;
        tm_array* log_likelihood;
        tm_array* log_posterior;
        tm_array* steps;
        tm_position* position;
        tm_position* next;
        double beta;
        rd_rng* rng;
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
        int type;
        int nchains;
        int nparameters;
        tm_chain** chains;
        tm_statistics* stats;
        tm_functions* functions;
        tm_options* options;
        void* user;
        int flags;
        int steps;
    } tm_state;

    struct _tm_prior;
    typedef struct _tm_prior tm_prior;

    typedef double (*tm_likelihood)(tm_array*, tm_state*);
    /*< The likelihood function. */
    typedef double (*tm_prior_function)(double, tm_prior*, tm_state*);
    /*< A prior function. */
    typedef void (*tm_proposal)(tm_position*, tm_position*, tm_state*, tm_chain*);
    /*< A proposal function. */
    typedef int (*tm_callback)(int, int, tm_state*);

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
        tm_likelihood log_likelihood;
        /**< The likelihood function (required) */
        tm_prior** priors;
        /**< Prior functions, specified for each parameter (required) */
        tm_proposal proposal;
        /**< Proposal function */
        tm_callback* callback;
        /**< A callback. */
        int nparameters;
    };

    struct _tm_options {
        int nchains;
        /**< Number of chains to use to check convergence on */
        double acceptance_rate_goal;
        /**< Goal acceptance rate (used to determine steps) */
        int discard;
        /**< Number of steps to discard to minimize correlations */
        int exchanges;
        /**< For PT, how often to try an exchange */
        uint32_t seed;
        /**< The initial seed for the simulation */
        void* user;
        /**< A user-defined pointer */
    };


    /* Prior functions */
    tm_prior* tm_prior_new(tm_prior_function function, double min, double max, void* user);
    tm_prior** tm_priors_new(int nparameters);

    tm_prior* tm_prior_uniform(double min, double max);
    tm_prior* tm_prior_modified_jeffreys(double max, double a);
    tm_prior* tm_prior_jeffreys(double min, double max);

    tm_functions* tm_functions_new(int nparameters);

    /* Options */
    tm_options* tm_options_new(int nparameters);


    /* Proposal functions */
    void tm_proposal_gaussian_gibbs1(tm_position* x0, tm_position* x1, tm_state* state, tm_chain* chain);

    /* Initialize algorithm */
    tm_state* tm_new(int type, int nparameters,
            tm_functions* functions,
            tm_array* steps,
            tm_options* options,
            tm_error* error);
    void tm_free(tm_state* state);

    /* Read/write state */
    tm_error tm_write(char* path_out, tm_state* state);
    tm_state* tm_read(char* path_in, tm_error* error);
    void tm_fprintf(FILE* out, tm_state* state);

    /* Actions */
    void tm_set_position(tm_state* state, int chain, double* x);
    void tm_set_seed(tm_state* state, int chain, uint32_t seed);

    void tm_determine_steps(tm_state* state);
    int tm_burn_in(tm_state* state, tm_error* error);
    int tm_run(tm_state* state, int steps, tm_error* error);

    /* Array allocation */
    tm_array* tm_array_new(double* array, const int len);
    void tm_array_free(tm_array* array);
    double tm_array_push(tm_array* array, const double v);
    tm_array* tm_array_clone(const tm_array* array);
    tm_array* tm_array_of(const int len, const double value);


#ifdef	__cplusplus
}
#endif

#endif	/* TINYMCMC_H */

