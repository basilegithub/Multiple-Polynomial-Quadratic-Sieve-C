#ifndef COMPUTE_SQRT_H
#define COMPUTE_SQRT_H

#include <gmp.h>
#include <stdbool.h>

#include "structures.h"

void build_sqrt(mpz_t n,unsigned long dim, bool vector[dim], dyn_array_classic primes, mpz_t x, mpz_t y, dyn_array relations, dyn_array smooth);

#endif // COMPUTE_SQRT_H