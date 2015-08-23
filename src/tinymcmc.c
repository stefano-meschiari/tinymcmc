#include "tinymcmc.h"
#include "stdlib.h"
#include "math.h"

#define TM_VALLOC(n) malloc(sizeof(double) * n)
#define TM_COPY(dest, src) do { for (int __i = 0; __i < dest->len; __i++) \
    dest->v[__i] = src->v[__i]; } while(0)

/**
 * Creates a custom prior function object.
 * @param function A pointer to a function 
 * @param min
 * @param max
 * @param user
 * @return 
 */
tm_prior* tm_new_prior(tm_prior_function function, double min, double max, void* user) {
    tm_prior* p = malloc(sizeof (tm_prior));
    p->min = min;
    p->max = max;
    p->prior_function = function;
    p->user = user;
    p->extra = NULL;
    return p;
}

/**
 * Frees the prior object.
 * @param prior Object to be freed
 */
void tm_free_prior(tm_prior* prior) {
    free(prior->extra);
    free(prior);
}

double _tm_uniform_f(double v, tm_prior* prior, tm_state* state) {
    return (v >= prior->min && v <= prior->max ? 1 : 0) / (prior->max - prior->min);
}

/**
 * A prior with uniform probability between min and max (inclusive).
 * @param min Minimum value 
 * @param max Maximum value
 * @return A prior object
 */
tm_prior* tm_uniform_prior(double min, double max) {
    return tm_new_prior(_tm_uniform_f, min, max, NULL);
}

double _tm_jeffreys_f(double v, tm_prior* prior, tm_state* state) {
    return 1. / ((v >= prior->min && v <= prior->max ? 1 : 0) * prior->extra[0]);
}

/**
 * A Jeffreys prior between min and max.
 * @param min Minimum value 
 * @param max Maximum value
 * @return A prior object
 */
tm_prior* tm_jeffreys_prior(double min, double max) {
    tm_prior* p = tm_new_prior(_tm_jeffreys_f, min, max, NULL);
    p->extra = TM_VALLOC(1);
    p->extra[0] = log(max / min);
    return p;
}

double _tm_modified_jeffreys_f(double v, tm_prior* prior, tm_state* state) {
    return 1. / ((v + prior->extra[1]) * prior->extra[0]);
}

/**
 * A modified Jeffreys prior that extends between 0 and max,
 * and is symmetric about 0. a is a scale parameter.
 * @param max Maximum value 
 * @param a Scale parameter
 * @return A prior object.
 */
tm_prior* tm_modified_jeffreys_prior(double max, double a) {
    tm_prior* p = tm_new_prior(_tm_jeffreys_f, 0, max, NULL);
    p->extra = TM_VALLOC(2);
    p->extra[0] = log((a + max) / a);
    p->extra[1] = a;
    return p;
}

/**
 * A Gibbs sampler. Returns a new position x1 based on the previous position x0
 * by randomly drawing an integer coordinate i and setting 
 * 
 * x1[i] = x0[i] + runif() * steps[i]
 * x1[j != i] = x0[j]
 * 
 * steps are specified per-chain during setup.
 * 
 * @param x0 Initial position
 * @param x1 Final position (out parameter)
 * @param state State of the MCMC algorithm
 * @param chain The chain we are proposing a new position for.
 */
void tm_gaussian_gibbs1_proposal(tm_array* x0, tm_array* x1, tm_state* state, tm_chain* chain) {
    int par = rd_rng_uintb(&chain->rng, x0->len);
    TM_COPY(x1, x0);
    x1->v[par] += rd_rng_exp(&chain->rng) * chain->steps->v[par];
}

double tm_log_merit(tm_state* state, tm_array* x, double beta) {
    double log_likelihood = state->functions->log_likelihood(x, state);
    double prior = 0;
    for (int i = 0; i < x->len; i++) {
        tm_prior* pf = state->functions->priors[i];
        prior += log10(pf->prior_function(x->v[i],
                pf, state));
    }
    return beta * log_likelihood + prior;
}

/**
 * Takes a single step by drawing a new proposal position, calculating the
 * merit (log L + log (sum_i priors_i)), and deciding whether to accept or
 * reject the new position.
 * 
 * @param state State of the MCMC algorithm.
 * @param nchain Number of chain we are stepping forward.
 * @return true if a new position has been accepted, false otherwise.
 */
bool tm_step(tm_state* state, int nchain, tm_error* error) {
    tm_chain* chain = state->chains[nchain];
    tm_gaussian_gibbs1_proposal(chain->position, chain->next, state, chain);
    double merit1 = tm_log_merit(state, chain->next, chain->beta);
    if (chain->n == 0)
        chain->merit = tm_log_merit(state, chain->position, chain->beta);
    double a = MIN(exp(merit1 - chain->merit), 1);
    double u = rd_rng_double(&chain->rng);

    ++chain->n;
    if (u < a && IS_FINITE(merit1)) {
        chain->merit = merit1;
        TM_COPY(chain->position, chain->next);
        return true;
    }

    if (!IS_FINITE(merit1) && error != NULL) {
        *error |= TM_INFINITY_ENCOUNTERED;
    }

    return false;
}

/**
 * 
 * @param state
 * @param steps
 * @param error
 * @return 
 */
int tm_run(tm_state* state, int steps, tm_error* error) {
    // TODO: Finish here
    int every = state->options->discard;

    for (int nchain = 0; nchain < state->nchains; nchain++) {
        for (int i = 0; i < steps; i++) {

        }
    }
    return 0;
}