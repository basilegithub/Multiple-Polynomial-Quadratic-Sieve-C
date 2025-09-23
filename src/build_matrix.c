#include <gmp.h>
#include <stdio.h>

#include "structures.h"

void build_sparse_matrix(dyn_array_classic* A, dyn_array relations, dyn_array_classic primes, unsigned long *nonzero, double *density, unsigned long *nb_lines, dyn_array_classic* rel_weight)
{
    for (unsigned long i = 0 ; i < relations.len ; i++) append_classic(rel_weight,0);
    *nonzero = 0;
    *nb_lines = 0;
    mpz_t tmp2;
    mpz_init(tmp2);
    for (unsigned long i = 0 ; i < relations.len ; i++)
    {
        if (mpz_cmp_ui(*(relations.start+i),0) < 0)
        {
            append_classic(A,i);
            *nonzero += 1;
            *(rel_weight->start+i) += 1;
        }
    }
    append_classic(A,relations.len);
    *nb_lines += 1;
    for (unsigned long i = 0 ; i < primes.len ; i++)
    {
        for (unsigned long j = 0 ; j < relations.len ; j++)
        {
            unsigned long tmp = 0;
            mpz_set_ui(tmp2,*(primes.start+i));
            while (mpz_divisible_p(*(relations.start+j),tmp2))
            {
                tmp ^= 1;
                mpz_mul_ui(tmp2,tmp2,*(primes.start+i));
            }
            if (tmp)
            {
                append_classic(A,j);
                *nonzero += 1;
                *(rel_weight->start+j) += 1;
            }
        }
        if (*(A->start+A->len-1) != relations.len)
        {
            append_classic(A,relations.len);
            *nb_lines += 1;
        }
    }
    *density = (double) (*nonzero)/(*nb_lines);
    mpz_clear(tmp2);
}

void build_dense_matrix(dyn_array relations, dyn_array_classic primes, unsigned long relations_len, unsigned long base_size, mpz_t *dense_matrix)
{

    mpz_t tmp, tmp2;
    mpz_inits(tmp, tmp2, NULL);

    unsigned long tmp_char;

    for (unsigned long i = 0 ; i < relations_len ; i++)
    {
        if (mpz_cmp_ui(*(relations.start+i), 0) < 0) mpz_set_ui(tmp, 1);
        else mpz_set_ui(tmp, 0);

        for (unsigned long j = 0 ; j < primes.len ; j++)
        {
            mpz_mul_2exp(tmp, tmp, 1);

            tmp_char = 0;
            mpz_set_ui(tmp2, *(primes.start+j));

            while (mpz_divisible_p(*(relations.start+i), tmp2))
            {
                tmp_char ^= 1;
                mpz_mul_ui(tmp2, tmp2, *(primes.start+j));
            }

            mpz_add_ui(tmp, tmp, tmp_char);
        }

        mpz_set(dense_matrix[i], tmp);

    }
    mpz_clears(tmp, tmp2, NULL);
}