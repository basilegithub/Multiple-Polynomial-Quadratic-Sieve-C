#ifndef REDUCE_MATRIX_H
#define REDUCE_MATRIX_H

#include <gmp.h>

#include "structures.h"

void reduce_relations(dyn_array* relations, dyn_array* smooth, dyn_array_classic* primes, mpz_t n);
void reduce_matrix(dyn_array_classic* matrix, dyn_array* relations, dyn_array* smooth, unsigned long limit, mpz_t n, dyn_array_classic* rel_weight, unsigned long merge_bound);

#endif // REDUCE_MATRIX_H