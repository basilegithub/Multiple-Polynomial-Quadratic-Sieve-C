#ifndef BUILD_MATRIX_H
#define BUILD_MATRIX_H

#include <gmp.h>

#include "structures.h"

void build_matrix(dyn_array_classic* A, dyn_array relations, dyn_array_classic primes, unsigned long *nonzero, double *density, unsigned long *nb_lines, dyn_array_classic* rel_weight);

#endif // BUILD_MATRIX_H