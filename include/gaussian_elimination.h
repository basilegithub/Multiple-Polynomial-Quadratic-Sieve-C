#ifndef GAUSSIAN_ELIMINATION_H
#define GAUSSIAN_ELIMINATION_H

#include <stdbool.h>

#include "structures.h"

void gaussian_elimination(mpz_t *dense_matrix, mpz_t *res, unsigned long relations_len, unsigned long base_size);
bool row_is_zero(mpz_t *dense_matrix, size_t row_index);

#endif // GAUSSIAN_ELIMINATION_H