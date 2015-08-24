#include "tinymcmc.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include "stdio.h"

#define TM_VALLOC(n) malloc(sizeof(double) * n)
#define TM_COPY(dest, src) do { for (int __i = 0; __i < dest->len; __i++) \
    dest->v[__i] = src->v[__i]; } while(0)


#define TM_MIN_CAPACITY 16
#define tm_unrecoverable_error(error) _tm_unrecoverable_error(error, __FILE__, __LINE__)

void _tm_unrecoverable_error(const char* error, const char* file, const int line) {
    fprintf(stderr, "[%s:%d] %s\n",
            file, line, error);
    abort();
}

/**
 * Allocates a new expandable array.
 * @param array Array to use as the backing store, or NULL.
 * @param len The length of the backing store array.
 * @return A new array
 */
tm_array* tm_array_new(double* array, const int len) {
    tm_array* a = malloc(sizeof (tm_array));
    if (a == NULL) {
        tm_unrecoverable_error("Could not allocate memory for array.");
    }
    if (array != NULL) {
        a->v = array;
        a->len = a->capacity = len;
    } else {
        a->v = malloc(sizeof (double) * TM_MIN_CAPACITY);
        a->len = 0;
        a->capacity = TM_MIN_CAPACITY;
    }
    return a;
}

/**
 * Frees the expandable array.
 * @param array
 */
void tm_array_free(tm_array* array) {
    if (array != NULL) {
        free(array->v);
        free(array);
    }
}

/**
 * Appends a number to the end of the array
 * @param array Expandable array
 * @param v New value to push
 * @return Pushed value
 */
double tm_array_push(tm_array* array, double v) {
    if (array->len == array->capacity) {
        array->capacity = 2 * array->capacity;
        array->v = realloc(array->v, array->capacity);
        if (array->v == NULL)
            tm_unrecoverable_error("Could not expand array.");
    }

    array->v[array->len] = v;
    array->len++;
    return v;
}

/**
 * Clones an array
 * @param array 
 * @return A cloned array
 */
tm_array* tm_array_clone(const tm_array* array) {
    double* arr = TM_VALLOC(array->len);
    memcpy(arr, array->v, sizeof (double) * array->len);
    return tm_array_new(arr, array->len);
}

tm_array* tm_array_of(const int len, const double value) {
    double* arr = TM_VALLOC(len);
    for (int i = 0; i < len; i++)
        arr[i] = value;
    return tm_array_new(arr, len);
}

void tm_array_set(tm_array* dest, const double* source) {
    for (int i = 0; i < dest->len; i++)
        dest->v[i] = source[i];
}

/**
 * Creates a custom prior function object.
 * @param function A pointer to a function 
 * @param min
 * @param max
 * @param user
 * @return 
 */
tm_prior* tm_prior_new(tm_prior_function function, double min, double max, void* user) {
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
void tm_prior_free(tm_prior* prior) {
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
tm_prior* tm_prior_uniform(double min, double max) {
    return tm_prior_new(_tm_uniform_f, min, max, NULL);
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
tm_prior* tm_prior_jeffreys(double min, double max) {
    tm_prior* p = tm_prior_new(_tm_jeffreys_f, min, max, NULL);
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
tm_prior* tm_prior_modified_jeffreys(double max, double a) {
    tm_prior* p = tm_prior_new(_tm_jeffreys_f, 0, max, NULL);
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
void tm_proposal_gaussian_gibbs1(tm_position* x0, tm_position* x1, tm_state* state, tm_chain* chain) {
    int par = rd_rng_uintb(chain->rng, x0->x->len);
    TM_COPY(x1->x, x0->x);
    x1->x->v[par] += rd_rng_exp(chain->rng) * chain->steps->v[par];
}

/**
 * Calculates the posterior for the given position,
 * p(theta | D) = L(D | theta) x p(theta)
 * @param state Current state.
 * @param position Position (the fields log_likelihood and log_posterior
 * are filled by this function)
 * @param beta Temperature
 */
void tm_log_posterior(tm_state* state, tm_position* position, const double beta) {
    double loglik = state->functions->log_likelihood(position->x, state);
    double prior = 0;
    for (int i = 0; i < position->x->len; i++) {
        tm_prior* pf = state->functions->priors[i];
        prior += log10(pf->prior_function(position->x->v[i],
                pf, state));
    }
    position->log_likelihood = loglik;
    position->log_posterior = beta * loglik + prior;
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
    tm_proposal_gaussian_gibbs1(chain->position, chain->next, state, chain);

    double loglik;
    tm_log_posterior(state, chain->next, chain->beta);
    if (chain->length == 0) {
        tm_log_posterior(state, chain->position, chain->beta);
    }
    double a = MIN(exp(chain->next->log_posterior - chain->position->log_posterior), 1);
    double u = rd_rng_double(chain->rng);

    ++chain->length;
    if (u < a && IS_FINITE(chain->next->log_posterior)) {
        chain->position->log_posterior = chain->next->log_posterior;
        chain->position->log_likelihood = chain->next->log_likelihood;

        TM_COPY(chain->position->x, chain->next->x);
        return true;
    }

    if (!IS_FINITE(chain->next->log_posterior) && error != NULL) {
        *error |= TM_INFINITY_ENCOUNTERED;
    }

    return false;
}

void tm_append_chain(tm_chain* chain) {
    for (int i = 0; i < chain->position->x->len; i++) {
        tm_array_push(chain->values[i], chain->position->x->v[i]);
    }
    tm_array_push(chain->log_likelihood, chain->position->log_likelihood);
    tm_array_push(chain->log_posterior, chain->position->log_posterior);
}

/**
 * Runs the MCMC chain for the specified number of steps.
 * @param state The state of the MCMC algorithm.
 * @param steps The number of steps.
 * @param error A possible error code.
 * @return The number of steps actually taken.
 */
int tm_run(tm_state* state, int steps, tm_error* error) {
    int every = state->options->discard;
    int exchanges = state->options->exchanges;

    int par_steps = (exchanges <= 0 ? steps : exchanges);
    int step = 0;

    while (step <= steps) {
        #pragma parallel for
        for (int nchain = 0; nchain < state->nchains; nchain++) {
            for (int i = 0; i < par_steps; i++) {
                tm_step(state, nchain, error);

                if (i + step % every == 0) {
                    tm_append_chain(state->chains[nchain]);
                }
            }
        }
        step += par_steps;
    }
    return steps;
}

/**
 * Creates a new MCMC state object.
 * @param type The type of algorithm (TM_BASE or TM_PARALLEL_TEMPERING)
 * @param nparameters Number of parameters
 * @param functions A pointer to a functions struct, containing the log-likelihood
 * function and the priors.
 * @param steps The initial steps for each parameter (will be optimized at the start of the
 * run)
 * @param options A pointer to the options struct.
 * @param error 
 * @return 
 */
tm_state* tm_new(int type, int nparameters,
        tm_functions* functions,
        tm_array* steps,
        tm_options* options,
        tm_error* error) {
    tm_state* s = malloc(sizeof (tm_state));
    s->functions = functions;
    s->options = options;
    s->type = type;
    s->nparameters = nparameters;
    s->flags = 0;

    s->user = s->options->user;
    s->nchains = options->nchains;
    s->stats = malloc(sizeof (tm_statistics));
    s->chains = malloc(s->nchains * sizeof (tm_chain*));

    for (int i = 0; i < s->nchains; i++) {
        s->chains[i] = malloc(sizeof (tm_chain));
        s->chains[i]->beta = 0;
        s->chains[i]->length = 0;

        s->chains[i]->position = malloc(sizeof (tm_position));
        s->chains[i]->position->x = tm_array_of(nparameters, 0.);

        s->chains[i]->next = malloc(sizeof (tm_position));
        s->chains[i]->next->x = tm_array_of(nparameters, 0.);

        s->chains[i]->rng = malloc(sizeof (rd_rng));
        rd_rng_init_time(s->chains[i]->rng);

        s->chains[i]->steps = tm_array_clone(steps);
        s->chains[i]->log_likelihood = tm_array_new(NULL, 0);
        s->chains[i]->log_posterior = tm_array_new(NULL, 0);

        s->chains[i]->values = malloc(nparameters * sizeof (tm_array));
    }

    return s;
}

/**
 * Sets the initial or current position for a chain.
 * @param state State of the MCMC algorithm
 * @param chain Chain ID to be modified
 * @param x Position to use (an array with at least nparameters elements)
 */
void tm_set_position(tm_state* state, int chain, double* x) {
    tm_array_set(state->chains[chain]->position->x, x);
}

tm_options* tm_options_new(int nparameters) {
    tm_options* options = malloc(sizeof (tm_options));
    options->acceptance_rate_goal = (nparameters > 1 ? 0.25 : 0.44);
    options->burn_in = nparameters * 500;
    options->discard = nparameters * 10;
    options->seed = 1234;
    options->nchains = 4;
    options->exchanges = -1;
    return options;
}

void tm_options_free(tm_options* options) {
    free(options);
}
