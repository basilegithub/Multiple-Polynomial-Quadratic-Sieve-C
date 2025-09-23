#ifndef BUILD_MATRIX_H
#define BUILD_MATRIX_H

#include "structures.h"

void build_sparse_matrix(dyn_array_classic* A, dyn_array relations, dyn_array_classic primes, unsigned long *nonzero, double *density, unsigned long *nb_lines, dyn_array_classic* rel_weight);
void build_dense_matrix(dyn_array relations, dyn_array_classic primes, unsigned long relations_len, unsigned long base_size, mpz_t *dense_matrix);

#endif // BUILD_MATRIX_H