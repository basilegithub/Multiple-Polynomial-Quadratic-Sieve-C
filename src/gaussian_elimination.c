#include <gmp.h>
#include <stdio.h>

void gaussian_elimination(unsigned long relations_len, unsigned long base_size, mpz_t *dense_matrix, mpz_t *res)
{
    mpz_t tmp, tmp2;
    mpz_inits(tmp, tmp2, NULL);

    unsigned long tmp_long;
    unsigned long i;
    int flag;

    mpz_t * restrict DM = dense_matrix;
    mpz_t * restrict R = res;

    mpz_set_ui(tmp, 1);

    for (unsigned long i = 0 ; i < relations_len ; i++)
    {
        mpz_mul_2exp(tmp2, tmp, relations_len-i-1);
        mpz_set(R[i], tmp2);
    }

    unsigned long index = 0;

    for (unsigned long j = 0 ; j < base_size ; j++)
    {
        i = index;
        
        mpz_div_2exp(tmp, DM[i], j);
        mpz_set_ui(tmp2, 1);
        mpz_and(tmp, tmp, tmp2);
        flag = (mpz_cmp_ui(tmp, 0) == 0);

        while (i < relations_len && flag)
        {
            i++;

            mpz_div_2exp(tmp, DM[i], j);
            mpz_set_ui(tmp2, 1);
            mpz_and(tmp, tmp, tmp2);
            flag = (mpz_cmp_ui(tmp, 0) == 0);
        }

        if (i < relations_len && !flag) // We have found the pivot row
        {
            if (i != index) // If needed, permute the rows so that pivot is at the top
            {
                mpz_set(tmp, DM[i]);
                mpz_set(DM[i], DM[index]);
                mpz_set(DM[index], tmp);

                mpz_set(tmp, R[i]);
                mpz_set(R[i], R[index]);
                mpz_set(R[index], tmp);
            }

            for (unsigned long k = i+1 ; k < relations_len ; k++) // reduce all the lower rows so that coefficient dense_matrix[k][j] = 0
            {
                mpz_div_2exp(tmp, DM[k], j);
                mpz_set_ui(tmp2, 1);
                mpz_and(tmp, tmp, tmp2);
                flag = (mpz_cmp_ui(tmp, 1) == 0);

                if (flag)
                {
                    mpz_xor(DM[k], DM[k], DM[index]);
                    mpz_xor(R[k], R[k], R[index]);
                }
            }
            index++;
        }
    }
    mpz_clears(tmp, tmp2, NULL);
}

unsigned char row_is_zero(unsigned long relations_len, unsigned long base_size, mpz_t *dense_matrix, unsigned long row_index)
{
    return (mpz_cmp_ui(dense_matrix[row_index], 0) == 0);
}