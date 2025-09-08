#ifndef GENERATE_PRIMES_H
#define GENERATE_PRIMES_H

#include <gmp.h>

#include "structures.h"

void smoothB(mpz_t B, dyn_array_classic* primes);

#endif // GENERATE_PRIMES_H