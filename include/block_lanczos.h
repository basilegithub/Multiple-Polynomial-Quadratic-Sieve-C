#ifndef BLOCK_LANCZOS_H
#define BLOCK_LANCZOS_H

#include "structures.h"

unsigned int switch_indices(size_t d, size_t mask);
void multiply_d(size_t *output, const size_t *dense_matrix, const size_t d, const size_t N);
void multiply_d_inplace(size_t *dense_matrix, const size_t d, const size_t N);
void extract_columns(size_t *W_inv, size_t *d, size_t *T, size_t N);
void solve(mpz_t *matrix, mpz_t *kernel, size_t nb_rows, size_t matrix_len);
void block_lanczos(dyn_array *output, dyn_array_classic sparse_matrix, size_t nb_relations, size_t block_size, unsigned long index, FILE *logfile);

#endif // BLOCK_LANCZOS_H