#include <gmp.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "structures.h"
#include "utils.h"

void poly_anul(mpz_t D, mpz_t B, unsigned long m)
{
    mpz_t A, C, E, Q, R, degB;

    mpz_init_set_ui(A, 1);
    mpz_mul_2exp(A, A, m<<1);
    mpz_init_set_ui(C, 0);
    mpz_set_ui(D, 1);

    mpz_inits(degB, E, Q, R, NULL);

    mpz_set(degB, B);
    my_int_log2(degB);

    while (mpz_cmp_ui(degB, m) > -1)
    {
        div_poly(Q, R, A, B);
        poly_prod(degB, Q, D);
        mpz_xor(E, C, degB);

        mpz_set(C, D);
        mpz_set(D, E);
        mpz_set(A, B);
        mpz_set(B, R);
        mpz_set(degB, B);
        my_int_log2(degB);
    }
    mpz_clears(A, C, E, Q, R, degB, NULL);
}

void wiedemann(dyn_array_classic A, mpz_t poly_res, unsigned long n, bool vec[n], unsigned long limit, unsigned long degree)
{
    bool vector_tmp[n];
    bool block[n];

    memset(vector_tmp, 0, n * sizeof *vector_tmp);
    memset(block, 0, n * sizeof *block);
    
    multiply(A, n, limit, vec, block);
    mpz_t poly_product, sequence, annihilator_poly;
    mpz_init_set_ui(poly_product, 1);
    mpz_inits(sequence, annihilator_poly, NULL);

    unsigned long d = degree;

    bool lambda[n];
    bool flag = true;
    while (flag)
    {
        for (size_t i = 0 ; i < n ; i++)
        {
            lambda[i] = rand()&1;
        }
        memcpy(vector_tmp, block, n * sizeof(bool));

        mpz_set_ui(sequence, 0);
        for (unsigned long i = 0 ; i < ((n-d)<<1) - 1 ; i++)
        {
            mpz_add_ui(sequence,sequence, dot_prod(n, lambda, vector_tmp));
            mpz_mul_2exp(sequence, sequence, 1);
            multiply(A, n, limit, vector_tmp, vector_tmp);
        }
        mpz_add_ui(sequence, sequence, dot_prod(n, lambda, vector_tmp));
        mpz_set_ui(annihilator_poly, 1);

        poly_anul(annihilator_poly, sequence, n-d);
        poly_prod(poly_product, poly_product, annihilator_poly);
        poly_eval(A, annihilator_poly, n, block, block, limit);
        my_int_log2(annihilator_poly);
        d += mpz_get_ui(annihilator_poly);

        flag = false;
        for (size_t i = 0 ; i < n && !flag; i++)
        {
            if (block[i]) flag = true;
        }
    }
    mpz_set(poly_res, poly_product);
    mpz_clears(poly_product, sequence, annihilator_poly, NULL);
}