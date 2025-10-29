#include <gmp.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "structures.h"
#include "utils.h"

void poly_anul(mpz_t D, mpz_t B, unsigned long m)
{
    mpz_t A, C, E, Q, R, tmp;

    mpz_init_set_ui(A, 1);
    mpz_mul_2exp(A, A, m<<1);
    mpz_init_set_ui(C, 0);
    mpz_set_ui(D, 1);

    mpz_inits(tmp, E, Q, R, NULL);

    while (mpz_sizeinbase(B, 2) > m)
    {
        div_poly(Q, R, A, B);
        poly_prod(tmp, Q, D);
        mpz_xor(E, C, tmp);

        mpz_set(C, D);
        mpz_set(D, E);
        mpz_set(A, B);
        mpz_set(B, R);
    }
    mpz_clears(A, C, E, Q, R, tmp, NULL);
}

void find_kernel_vectors(dyn_array *kernel_vectors, dyn_array_classic A, mpz_t minimal_polynomial_estimate, size_t *block, size_t block_size, unsigned long n, unsigned long limit)
{
    mpz_t tmp_poly;
    mpz_init_set(tmp_poly, minimal_polynomial_estimate);

    size_t cpt = 2;
    while (mpz_even_p(tmp_poly))
    {
        cpt++;
        mpz_div_2exp(tmp_poly, tmp_poly, 1);
    }

    size_t *res = calloc(n ,sizeof(size_t));
    // memcpy(res, block, n*sizeof(size_t));

    size_t *tmp_array = calloc(n, sizeof(size_t));

    size_t deg = mpz_sizeinbase(tmp_poly, 2)-1;

    for (size_t i = 0 ; i < deg ; i++)
    {
        if (mpz_tstbit(tmp_poly, deg-i))
        {
            for (size_t j = 0 ; j < n ; j++)
            {
                res[j] ^= block[j];
            }
        }
        multiply_size_t(A, n, limit, res, tmp_array);
        memcpy(res, tmp_array, n*sizeof(size_t));
    }

    if (mpz_tstbit(tmp_poly, 0))
    {
        for (size_t j = 0 ; j < n ; j++)
        {
            res[j] ^= block[j];
        }
    }

    // size_t *image = calloc(n, sizeof(size_t));
    // multiply_size_t(A, n, limit, res, image); // image = A.res

    size_t *indexes = calloc(block_size, sizeof(size_t));
    for (size_t i = 0 ; i < block_size ; i++)
    {
        indexes[i] = i;
    }

    unsigned int indexes_size = block_size;

    mpz_t kernel_vec;
    mpz_init(kernel_vec);

    for (size_t i = 0 ; i < cpt ; i++)
    {
        multiply_size_t(A, n, limit, res, tmp_array);

        for (size_t j = 0 ; j < indexes_size ; j++)
        {
            bool flag_is_zero = true;
            for (size_t k = 0 ; k < n ; k++)
            {
                if ((tmp_array[k]>>(block_size-indexes[j]-1))&1)
                {
                    flag_is_zero = false;
                    break;
                }
            }

            if (flag_is_zero)
            {
                mpz_set_ui(kernel_vec, 0);
                for (size_t k = 0 ; k < n-1 ; k++)
                {
                    mpz_add_ui(kernel_vec, kernel_vec, (res[k]>>(block_size-indexes[j]-1))&1);
                    mpz_mul_2exp(kernel_vec, kernel_vec, 1);
                }
                mpz_add_ui(kernel_vec, kernel_vec, (res[n-1]>>(block_size-indexes[j]-1))&1);
                append(kernel_vectors, kernel_vec);

                indexes[j] = indexes[--indexes_size];
            }
        }
        memcpy(res, tmp_array, n*sizeof(size_t));
    }

    mpz_clears(kernel_vec, tmp_poly, NULL);

    free(res);
    free(tmp_array);
    free(indexes);
}

void wiedemann(dyn_array *kernel_vectors, dyn_array_classic A, mpz_t minimal_polynomial_estimate, size_t block_size, unsigned long n, unsigned long limit)
{
    // create a new block
    size_t *block = calloc(n, sizeof(size_t));

    size_t *tmp_block = calloc(n, sizeof(size_t));
    size_t *tmp_block2 = calloc(n, sizeof(size_t));

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

    // multiply_size_t(A, n, limit, block, tmp_block);

    memcpy(tmp_block, block, n*sizeof(size_t));

    mpz_t *sequences = malloc(block_size*sizeof(mpz_t));

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
            size_t shifted = (dot_prod_value >> (block_size - j - 1));
            size_t bit_result;
            
            // Force the AND operation via assembly
            __asm__ ("and %1, %0" : "=r" (bit_result) : "r" (shifted), "0" (1));
            // For some unkown reason it does not work in pure C
            // It just does bit_result = shifted&1 ;
            
            mpz_add_ui(sequences[j], sequences[j], bit_result);
            mpz_mul_2exp(sequences[j], sequences[j], 1);
        }

        multiply_size_t(A, n, limit, tmp_block, tmp_block2);
        memcpy(tmp_block, tmp_block2, n*sizeof(size_t));
    }

    dot_prod_value = dot_prod(n, lambda, tmp_block);
    for (size_t j = 0 ; j < block_size ; j++)
    {
        mpz_add_ui(sequences[j], sequences[j], (dot_prod_value>>(block_size-j-1))&1);
    }

    // compute the annealing polynomial of each sequence

    mpz_t anneal_polynomial, gcd, q, r;
    mpz_inits(anneal_polynomial, gcd, q, r, NULL);

    size_t deg;

    for (size_t i = 0 ; i < block_size ; i++)
    {
        mpz_set_ui(anneal_polynomial, 1);
        mpz_set(r, sequences[i]);
        poly_anul(anneal_polynomial, sequences[i], n);
        deg = mpz_sizeinbase(anneal_polynomial, 2)-1;
        

        if (deg) // update the minimal polynomial
        {
            // compute gcd
            gcd_poly(gcd, anneal_polynomial, minimal_polynomial_estimate);
            // divide new poly by gcd
            div_poly(q, r, anneal_polynomial, gcd);
            // multiply estimate by the result
            poly_prod(minimal_polynomial_estimate, minimal_polynomial_estimate, q);
        }
    }


    mpz_clears(anneal_polynomial, gcd, q, r, NULL);

    // use new estimate to find kernel vectors
    find_kernel_vectors(kernel_vectors, A, minimal_polynomial_estimate, block, block_size, n, limit);

    free(block);
    free(tmp_block);
    free(tmp_block2);
    free(lambda);

    for (size_t i = 0 ; i < block_size ; i++) mpz_clear(sequences[i]);
    free(sequences);
}