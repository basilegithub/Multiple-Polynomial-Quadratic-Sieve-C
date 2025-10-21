#ifndef REDUCE_MATRIX_H
#define REDUCE_MATRIX_H

#include <gmp.h>

#include "structures.h"

void reduce_relations(dyn_array* relations, dyn_array* smooth, dyn_array_classic* primes, mpz_t n);
void bubble_sort_down(dyn_array_classic *sorted, dyn_array_classic weights, unsigned long tmp3);
void bubble_sort_up(dyn_array_classic *sorted, dyn_array_classic weights, unsigned long tmp3);
size_t get_index(dyn_array_classic sorted, unsigned long line_index);
bool delete_empty_row(dyn_array_classic *matrix, dyn_array_classic *weights, dyn_array_classic *sorted, size_t i, unsigned long row_index);

void delete_singleton(
    dyn_array *relations,
    dyn_array *smooth,
    dyn_array_classic *matrix,
    dyn_array_classic *weights,
    dyn_array_classic *rel_weight,
    dyn_array_classic sorted,
    size_t i,
    unsigned long row_index,
    unsigned long row_delimiter
);

void delete_two_elements_row(
    dyn_array *relations,
    dyn_array *smooth,
    dyn_array_classic *matrix,
    dyn_array_classic *rel_weight,
    dyn_array_classic *weights,
    dyn_array_classic sorted,
    mpz_t n,
    size_t i,
    unsigned long row_index,
    unsigned long row_delimiter
);

void reduce_matrix(
    dyn_array* relations,
    dyn_array* smooth,
    dyn_array_classic* matrix,
    dyn_array_classic* rel_weight,
    mpz_t n,
    unsigned long row_delimiter,
    unsigned long merge_bound
);

#endif // REDUCE_MATRIX_H