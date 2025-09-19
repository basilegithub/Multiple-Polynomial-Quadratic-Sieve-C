#ifndef SMOOTH_TEST_H
#define SMOOTH_TEST_H

#include <gmp.h>

#include "structures.h"

void batch_smooth(dyn_array* reported, dyn_array* large_primes, dyn_array_classic* smooth, dyn_array* tmp_array, mpz_t limit, mpz_t prod_primes, unsigned long prime);
void naive_smooth(dyn_array* reported, dyn_array* large_primes, dyn_array_classic* smooth, dyn_array_classic primes, mpz_t limit);

#endif // SMOOTH_TEST_H