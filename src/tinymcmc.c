#include "tinymcmc.h"
#include "stdlib.h"
#include "math.h"

#define tm_valloc(n) malloc(sizeof(double) * n)

tm_prior* tm_new_prior(tm_function function, double min, double max, void* user) {
    tm_prior* p = malloc(sizeof (tm_prior));
    p->min = min;
    p->max = max;
    p->prior_function = function;
    p->user = user;
    p->extra = NULL;
}

void tm_free_prior(tm_prior* prior) {
    free(prior->extra);
    free(prior);
}

double _tm_uniform_f(double v, tm_prior* prior, tm_state* state) {
    return (v >= prior->min && v <= prior->max ? 1 : 0) / (prior->max - prior->min);
}

tm_prior* tm_uniform_prior(double min, double max) {
    return tm_new_prior(_tm_uniform_f, min, max, NULL);
}

double _tm_jeffreys_f(double v, tm_prior* prior, tm_state* state) {
    return 1. / ((v >= prior->min && v <= prior->max ? 1 : 0) * prior->extra[0]);
}

tm_prior* tm_jeffreys_prior(double min, double max) {
    tm_prior* p = tm_new_prior(_tm_jeffreys_f, min, max, NULL);
    p->extra = tm_valloc(1);
    p->extra[0] = log(max / min);
    return p;
}

double _tm_modified_jeffreys_f(double v, tm_prior* prior, tm_state* state) {
    return 1. / ((v + prior->extra[1]) * prior->extra[0]);
}

tm_prior* tm_modified_jeffreys_prior(double max, double a) {
    tm_prior* p = tm_new_prior(_tm_jeffreys_f, 0, max, NULL);
    p->extra = tm_valloc(2);
    p->extra[0] = log((a + max) / a);
    p->extra[1] = a;
    return p;
}

