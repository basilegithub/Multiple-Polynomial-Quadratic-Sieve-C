#ifndef WIEDEMANN_H
#define WIEDEMANN_H

#include <gmp.h>
#include <stdbool.h>

#include "structures.h"

void poly_anul(mpz_t D, mpz_t B, unsigned long m);
void wiedemann(dyn_array *kernel_vectors, dyn_array_classic A, mpz_t minimal_polynomial_estimate, size_t block_size, unsigned long n, unsigned long limit);

#endif // WIEDEMANN_H