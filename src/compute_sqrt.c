#include <gmp.h>
#include <stdlib.h>
#include <stdbool.h>

#include "structures.h"

void build_sqrt(mpz_t n, unsigned long dim, bool vector[dim], dyn_array_classic factor_base, mpz_t x, mpz_t y, dyn_array relations, dyn_array smooth)
{
    mpz_t tmp_mpz, tmp_mpz_2;
    mpz_inits(tmp_mpz, tmp_mpz_2, NULL);

    mpz_set_ui(y, 1);
    mpz_set_ui(x, 1);
    unsigned long *vector_y = calloc(factor_base.len, sizeof(unsigned long));

    for (size_t i = 0 ; i < factor_base.len ; i++) vector_y[i] = 0;

    for (size_t e = 0 ; e < dim ; e++)
    {
        if (vector[e])
        {
            mpz_mul(tmp_mpz, x, smooth.start[e]);
            mpz_mod(x, tmp_mpz, n);

            for (size_t i = 0 ; i < factor_base.len ; i++)
            {
                mpz_set_ui(tmp_mpz, factor_base.start[i]);
                vector_y[i] += mpz_remove(tmp_mpz_2, relations.start[e], tmp_mpz);
            }
        }
    }
    
    for (size_t k = 0 ; k < factor_base.len ; k++)
    {
        mpz_set_ui(tmp_mpz, factor_base.start[k]);
        mpz_powm_ui(tmp_mpz, tmp_mpz, vector_y[k]>>1, n);
        mpz_mul(tmp_mpz_2, y, tmp_mpz);
        mpz_mod(y, tmp_mpz_2, n);
    }
    mpz_clears(tmp_mpz, tmp_mpz_2, NULL);
    free(vector_y);
}