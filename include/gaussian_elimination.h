#ifndef GAUSSIAN_ELIMINATION_H
#define GAUSSIAN_ELIMINATION_H

#include <stdbool.h>

#include "structures.h"

void gaussian_elimination(unsigned long relations_len, unsigned long base_size, mpz_t *dense_matrix, mpz_t *res);
bool row_is_zero(unsigned long relations_len, unsigned long base_size, mpz_t *dense_matrix, unsigned long row_index);

#endif // GAUSSIAN_ELIMINATION_H