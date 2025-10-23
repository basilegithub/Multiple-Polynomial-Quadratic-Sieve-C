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

void extract_columns(size_t *W_inv, size_t *d, size_t *T, size_t N)
{
    size_t *M = calloc(N<<1, sizeof(size_t));

    size_t *ident = calloc(N, sizeof(size_t));

    identity(ident, N);
    concatenate(M, T, ident, N);

    size_t *S = calloc(N, sizeof(size_t));
    size_t s_len = 0;

    for (size_t j = 0 ; j < N ; j++)
    {
        for (size_t k = j ; k < N ; k++)
        {
            if ((M[k]>>(2*N-j-1))&1)
            {
                size_t tmp = M[k];
                M[k] = M[j];
                M[j] = tmp;
                break;
            }
        }

        if ((M[j]>>(2*N-j-1))&1)
        {
            for (size_t k = 0 ; k < N ; k++)
            {
                if (k != j)
                {
                    M[k] ^= M[j];
                }
            }
            S[s_len++] = j;
        }
        else
        {
            for (size_t k = j ; k < N ; k++)
            {
                if ((M[k]>>(N-j-1))&1)
                {
                    size_t tmp = M[k];
                    M[k] = M[j];
                    M[j] = tmp;
                    break;
                }
            }

            if (!(M[j]>>(N-j-1))&1)
            {
                printf("you are fucked man\n\n\n");
            }

            for (size_t k = 0 ; k < N ; k++)
            {
                if (k != j)
                {
                    if ((M[k]>>(N-j-1))&1)
                    {
                        M[k] ^= M[j];
                    }
                }
            }
            M[j] = 0;
        }
    }

    for (size_t i = 0 ; i < N ; i++)
    {
        W_inv[i] = M[i]&((1<<N)-1);
    }

    *d = 0;
    for (size_t i = 0 ; i < s_len ; i++)
    {
        *d |= (1<<S[i]);
    }
}

void block_lanczos(dyn_array_classic sparse_matrix, size_t nb_relations, size_t block_size)
{
    size_t *Y = calloc(nb_relations, sizeof(size_t));

    for (size_t i = 0 ; i < nb_relations ; i++)
    {
        Y[i] = rand()%((1<<block_size)-1);
    }
}

void solve(mpz_t *matrix, mpz_t *kernel, size_t nb_relations, size_t matrix_len)
{
    size_t k = 0;

    for (size_t i = 0 ; i < nb_relations ; i++)
    {
        size_t index;
        for (size_t j = k ; j < matrix_len ; j++)
        {
            if (mpz_tstbit(matrix[j], nb_relations-i-1))
            {
                mpz_swap(matrix[j], matrix[k]);
                mpz_swap(kernel[j], kernel[k]);
                k++;
                break;
                index = j;
            }
        }

        for (size_t j = index+1 ; j < matrix_len ; j++)
        {
            if (mpz_tstbit(matrix[j], nb_relations-i-1))
            {
                mpz_xor(matrix[j], matrix[j], matrix[k-1]);
                mpz_xor(kernel[j], kernel[j], kernel[k-1]);
            }
        }
    }
}