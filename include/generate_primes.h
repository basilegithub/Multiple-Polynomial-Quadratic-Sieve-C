#ifndef GENERATE_PRIMES_H
#define GENERATE_PRIMES_H

#include <gmp.h>

#include "structures.h"

void erasthotenes_sieve(dyn_array_classic* primes, mpz_t bound);

#endif // GENERATE_PRIMES_H