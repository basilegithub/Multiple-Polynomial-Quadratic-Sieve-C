#include <gmp.h>
#include <stdio.h>

#include "structures.h"

void build_sparse_matrix(dyn_array relations, dyn_array_classic* matrix, dyn_array_classic* rel_weight, dyn_array_classic factor_base, unsigned long *nonzero_count, unsigned long *row_count, double *density)
{
    mpz_t current_prime, tmp_mpz;
    mpz_inits(current_prime, tmp_mpz, NULL);

    *nonzero_count = 0;
    *row_count = 0;

    for (size_t i = 0 ; i < relations.len ; i++) append_classic(rel_weight, 0);

    for (size_t i = 0 ; i < relations.len ; i++)
    {
        if (mpz_cmp_ui(relations.start[i], 0) < 0)
        {
            append_classic(matrix, i);
            (*nonzero_count)++;
            rel_weight->start[i]++;
        }
    }

    append_classic(matrix, relations.len);
    (*row_count)++;

    bool parity_bit;
    unsigned int bit_count;

    for (size_t i = 0 ; i < factor_base.len ; i++)
    {
        mpz_set_ui(current_prime, factor_base.start[i]);

        for (size_t j = 0 ; j < relations.len ; j++)
        {
            bit_count = mpz_remove(tmp_mpz, relations.start[j], current_prime);
            parity_bit = bit_count&1;

            if (parity_bit)
            {
                append_classic(matrix, j);
                (*nonzero_count)++;
                rel_weight->start[j]++;
            }
        }
        if (*(matrix->start + matrix->len - 1) != relations.len)
        {
            append_classic(matrix, relations.len);
            (*row_count)++;
        }
    }
    *density = (double) (*nonzero_count)/(*row_count);
    mpz_clears(current_prime, tmp_mpz, NULL);
}

void build_dense_matrix(dyn_array relations, dyn_array_classic factor_base, mpz_t *dense_matrix, unsigned long relations_len)
{

    mpz_t accumulator, current_prime, tmp_mpz;
    mpz_inits(accumulator, current_prime, tmp_mpz, NULL);

    mpz_t * restrict DM = dense_matrix;
    mpz_t * restrict RELS = relations.start;
    unsigned long * restrict PRIMES = factor_base.start;

    bool parity_bit;
    unsigned int bit_count;

    for (size_t i = 0 ; i < relations_len ; i++)
    {
        if (mpz_cmp_ui(RELS[i], 0) < 0) mpz_set_ui(accumulator, 1);
        else mpz_set_ui(accumulator, 0);

        for (size_t j = 0 ; j < factor_base.len ; j++)
        {
            mpz_set_ui(current_prime, PRIMES[j]);
            mpz_mul_2exp(accumulator, accumulator, 1);

            bit_count = mpz_remove(tmp_mpz, RELS[i], current_prime);
            parity_bit = bit_count&1;

            mpz_add_ui(accumulator, accumulator, parity_bit);
        }

        mpz_set(DM[i], accumulator);

    }
    mpz_clears(accumulator, current_prime,  tmp_mpz, NULL);
}