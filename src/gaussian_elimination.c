#include <gmp.h>
#include <stdio.h>
#include <stdbool.h>

void gaussian_elimination(mpz_t *dense_matrix, mpz_t *res, unsigned long relations_len, unsigned long base_size)
{
    mpz_t tmp;
    mpz_init(tmp);

    unsigned long row_index;
    int is_zero_bit;

    mpz_t * restrict matrix = dense_matrix;
    mpz_t * restrict result = res;

    mpz_set_ui(tmp, 1);

    for (size_t i = 0 ; i < relations_len ; i++)
    {
        mpz_mul_2exp(result[i], tmp, relations_len-i-1);
    }

    unsigned long pivot_row = 0;

    for (size_t j = 0 ; j < base_size ; j++)
    {
        row_index = pivot_row;
        is_zero_bit = (!mpz_tstbit(matrix[row_index], j));

        while (row_index < relations_len && is_zero_bit)
        {
            row_index++;
            is_zero_bit = (!mpz_tstbit(matrix[row_index], j));
        }

        if (row_index < relations_len && !is_zero_bit) // We have found the pivot row
        {
            if (row_index != pivot_row) // If needed, permute the rows so that pivot is at the top
            {
                mpz_swap(matrix[row_index], matrix[pivot_row]);
                mpz_swap(result[row_index], result[pivot_row]);
            }

            for (size_t k = row_index+1 ; k < relations_len ; k++) // reduce all the lower rows so that coefficient dense_matrix[k][j] = 0
            {
                is_zero_bit = (!mpz_tstbit(matrix[k], j));

                if (!is_zero_bit)
                {
                    mpz_xor(matrix[k], matrix[k], matrix[pivot_row]);
                    mpz_xor(result[k], result[k], result[pivot_row]);
                }
            }
            pivot_row++;
        }
    }
    mpz_clear(tmp);
}

bool row_is_zero(mpz_t *dense_matrix, size_t row_index)
{
    return ((bool) !mpz_cmp_ui(dense_matrix[row_index], 0));
}