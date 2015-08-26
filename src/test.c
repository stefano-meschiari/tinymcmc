#include "tinymcmc.h"
#include "math.h"

double gaussian_likelihood(tm_array* x, tm_state* state) {
    if (x->len == 1) {
        return -(x->v[0] * x->v[0]) / 2. - sqrt(2. * M_PI);
    } else {
        double r = 0;
        for (int i = 0; i < x->len; i++)
            r += -(x->v[i] * x->v[i]) / 2. - sqrt(2. * M_PI);
        return r;
    }
}

void test_gaussian(const int N) {

    tm_prior** priors = tm_priors_new(N);
    for (int i = 0; i < N; i++)
        priors[i] = tm_prior_uniform(-10, 10);

    tm_array* position = tm_array_of(N, 0.);

    tm_array* steps = tm_array_of(N, 0.1);

    tm_functions * f = tm_functions_new(N);
    f->log_likelihood = gaussian_likelihood;
    f->priors = priors;

    tm_options* options = tm_options_new(N);
    options->discard = 100;

    tm_error error;

    tm_state* s = tm_new(TM_MCMC_BASE,
                         N,
                         f,
                         steps,
                         options,
                         &error);
    for (int chain = 0; chain < s->nchains; chain++)
        tm_set_position(s, chain, position->v);

    tm_run(s, 1000000, &error);
    char name[512];
    sprintf(name, "test_gaussian_%d", N);

    tm_write(name, s);

    tm_free(s);
}

int main() {
    // Sample a one-dimensional gaussian
    printf("Running 1-d gaussian test...\n");
    test_gaussian(1);
    printf("Running 2-d gaussian test...\n");
    test_gaussian(2);
}
