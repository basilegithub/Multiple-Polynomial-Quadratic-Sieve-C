#include <gmp.h>
#include <stdio.h>
#include <stdbool.h>

void gaussian_elimination(unsigned long relations_len, unsigned long base_size, mpz_t *dense_matrix, mpz_t *res)
{
    mpz_t tmp, tmp2;
    mpz_inits(tmp, tmp2, NULL);

    unsigned long tmp_long;
    unsigned long i;
    int is_zero_bit;

    mpz_t * restrict matrix = dense_matrix;
    mpz_t * restrict result = res;

    mpz_set_ui(tmp, 1);

    for (size_t i = 0 ; i < relations_len ; i++)
    {
        mpz_mul_2exp(tmp2, tmp, relations_len-i-1);
        mpz_set(result[i], tmp2);
    }

    unsigned long pivot_row = 0;

    for (size_t j = 0 ; j < base_size ; j++)
    {
        i = pivot_row;
        
        mpz_div_2exp(tmp, matrix[i], j);
        mpz_set_ui(tmp2, 1);
        mpz_and(tmp, tmp, tmp2);
        is_zero_bit = (!mpz_cmp_ui(tmp, 0));

        mpz_set_ui(tmp2, 1);

        while (i < relations_len && is_zero_bit)
        {
            i++;

            mpz_div_2exp(tmp, matrix[i], j);
            mpz_and(tmp, tmp, tmp2);
            is_zero_bit = (!mpz_cmp_ui(tmp, 0));
        }

        if (i < relations_len && !is_zero_bit) // We have found the pivot row
        {
            if (i != pivot_row) // If needed, permute the rows so that pivot is at the top
            {
                mpz_set(tmp, matrix[i]);
                mpz_set(matrix[i], matrix[pivot_row]);
                mpz_set(matrix[pivot_row], tmp);

                mpz_set(tmp, result[i]);
                mpz_set(result[i], result[pivot_row]);
                mpz_set(result[pivot_row], tmp);
            }

            for (size_t k = i+1 ; k < relations_len ; k++) // reduce all the lower rows so that coefficient dense_matrix[k][j] = 0
            {
                mpz_div_2exp(tmp, matrix[k], j);
                mpz_set_ui(tmp2, 1);
                mpz_and(tmp, tmp, tmp2);
                is_zero_bit = (!mpz_cmp_ui(tmp, 1));

                if (is_zero_bit)
                {
                    mpz_xor(matrix[k], matrix[k], matrix[pivot_row]);
                    mpz_xor(result[k], result[k], result[pivot_row]);
                }
            }
            pivot_row++;
        }
    }
    mpz_clears(tmp, tmp2, NULL);
}

bool row_is_zero(unsigned long relations_len, unsigned long base_size, mpz_t *dense_matrix, size_t row_index)
{
    return ((bool) !mpz_cmp_ui(dense_matrix[row_index], 0));
}