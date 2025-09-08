#ifndef SMOOTH_TEST_H
#define SMOOTH_TEST_H

#include <gmp.h>

#include "structures.h"

void batch_smooth(dyn_array* reported, mpz_t prod_primes, mpz_t limit, dyn_array* large_primes, dyn_array_classic* smooth, dyn_array* tmp_array, unsigned long prime);

#endif // SMOOTH_TEST_H