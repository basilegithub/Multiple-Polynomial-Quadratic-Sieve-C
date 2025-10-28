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

void wiedemann(dyn_array *kernel_vectors, dyn_array_classic A, mpz_t poly_res, size_t block_size, unsigned long n, unsigned long limit, unsigned long degree)
{
    // create a new block
    size_t *block = calloc(n, sizeof(size_t));

    size_t tmp_block = calloc(n, sizeof(size_t));
    size_t tmp_block2 = calloc(n, sizeof(size_t));

    for (size_t i = 0 ; i < n ; i++)
    {
        for (size_t j = 0 ; j < block_size-1 ; j++)
        {
            block[i] ^= rand()&1;
            block[i] <<= 1;
        }
        block[i] ^= rand()&1;
    }

    // create one lbd projection

    bool *lambda = calloc(n, sizeof(bool));

    for (size_t i = 0 ; i < n ; i++)
    {
        lambda[i] = rand()&1;
    }

    // compute the sequences

    memcpy(tmp_block, block, n*sizeof(size_t));

    mpz_t *sequences = calloc(block_size, sizeof(mpz_t));
    for (size_t i = 0 ; i < block_size ; i++)
    {
        mpz_init_set_ui(sequences[i], 0);
    }

    size_t dot_prod_value;

    for (size_t i = 0 ; i < (2*n-1) ; i++)
    {
        dot_prod_value = dot_prod(n, lambda, tmp_block);
        for (size_t j = 0 ; j < block_size ; j++)
        {
            mpz_add_ui(sequences[j], sequences[j], (dot_prod_value>>(block_size-j-1))&1);
            mpz_mul_2exp(sequences[j], sequences[j], 1);
        }
        multiply(A, n, limit, tmp_block, tmp_block2);
        memcpy(tmp_block, tmp_block2, n*sizeof(size_t));
    }
    dot_prod_value = dot_prod(n, lambda, tmp_block);
    for (size_t j = 0 ; j < block_size ; j++)
    {
        mpz_add_ui(sequences[j], sequences[j], (dot_prod_value>>(block_size-j-1))&1);
    }

    // compute the annealing polynomial of each sequence

    mpz_t anneal_polynomial;
    mpz_init(anneal_polynomial);

    size_t deg;

    for (size_t i = 0 ; i < block_size ; i++)
    {
        poly_anul(anneal_polynomial, sequences[i], n);
        deg = mpz_sizeinbase(anneal_polynomial, 2)-1;

        if (deg) // update the minimal polynomial
        {

        }
    }

    mpz_clear(anneal_polynomial);

    // evaluate the whole block with the updated minimal polynomial
    // use the usual strat to find the kernel vectors
}