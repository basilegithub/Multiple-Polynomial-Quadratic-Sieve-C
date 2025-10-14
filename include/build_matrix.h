#ifndef BUILD_MATRIX_H
#define BUILD_MATRIX_H

#include "structures.h"

void build_sparse_matrix(dyn_array_classic* matrix, dyn_array relations, dyn_array_classic factor_base, unsigned long *nonzero_count, double *density, unsigned long *row_count, dyn_array_classic* rel_weight);
void build_dense_matrix(dyn_array relations, dyn_array_classic factor_base, unsigned long relations_len, unsigned long base_size, mpz_t *dense_matrix);

#endif // BUILD_MATRIX_H