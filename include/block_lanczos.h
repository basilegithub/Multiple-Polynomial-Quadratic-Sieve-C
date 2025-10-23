#ifndef BLOCK_LANCZOS_H
#define BLOCK_LANCZOS_H

#include "structures.h"

unsigned int switch_indices(size_t d, size_t mask);
void multiply_d(size_t *output, const size_t *dense_matrix, const size_t d, const size_t N);
void multiply_d_inplace(size_t *dense_matrix, const size_t d, const size_t N);

#endif // BLOCK_LANCZOS_H