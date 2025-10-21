#ifndef COMPUTE_SQRT_H
#define COMPUTE_SQRT_H

#include <gmp.h>
#include <stdbool.h>

#include "structures.h"

void build_sqrt(dyn_array relations, dyn_array smooth, dyn_array_classic factor_base, mpz_t n, mpz_t x, mpz_t y, unsigned long dim, bool vector[dim]);

#endif // COMPUTE_SQRT_H