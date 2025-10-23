#include <gmp.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "structures.h"
#include "utils.h"

unsigned int switch_indices(size_t d, size_t mask)
{
    return (~d) & mask;
}

void multiply_d(size_t *output, const size_t *dense_matrix, const size_t d, const size_t N)
{
    for (size_t i = 0 ; i < N ; i++)
    {
        output[i] = dense_matrix[i]&d;
    }
}

void multiply_d_inplace(size_t *dense_matrix, const size_t d, const size_t N)
{
    for (size_t i = 0 ; i < N ; i++)
    {
        dense_matrix[i] = dense_matrix[i]&d;
    }
}