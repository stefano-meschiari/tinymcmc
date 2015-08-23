#include <assert.h>
#include "../src/rng.c"

int main() {
  rd_rng r;
  rd_rng_init_time(&r);

  int i = 0;
  for (; i < 1024; ++i) {
    r.has_next = (i % 2 == 0);
    rd_rng_normal(&r);
  }
  return EXIT_SUCCESS;
}

