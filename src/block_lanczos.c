#include <gmp.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "structures.h"
#include "utils.h"

#define MIN(x, y) ((x) < (y) ? (x) : (y))

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

void block_lanczos(dyn_array *output, dyn_array_classic sparse_matrix, size_t nb_relations, size_t block_size, unsigned long index)
{
    size_t *Y = calloc(nb_relations, sizeof(size_t));

    for (size_t i = 0 ; i < nb_relations ; i++)
    {
        Y[i] = rand()%((1<<block_size)-1);
    }

    size_t *X = calloc(nb_relations, sizeof(size_t));

    size_t *tmp = calloc(nb_relations, sizeof(size_t));
    size_t *tmp2 = calloc(nb_relations, sizeof(size_t));

    size_t *V0 = calloc(nb_relations, sizeof(size_t));

    multiply_sparse(sparse_matrix, nb_relations, index, Y, tmp);
    sparse_multiply_transpose(sparse_matrix, tmp, V0, index);

    size_t *P = calloc(nb_relations, sizeof(size_t));

    size_t *V = calloc(nb_relations, sizeof(size_t));
    memcpy(V, V0, nb_relations*sizeof(size_t));

    size_t d = 1;
    size_t neg_d;
    size_t mask = (1<<block_size)-1;

    size_t i = 0;

    size_t *Z = calloc(nb_relations, sizeof(size_t));
    size_t *vAv = calloc(block_size, sizeof(size_t));
    size_t *vAAv = calloc(block_size, sizeof(size_t));

    size_t *W_inv = calloc(block_size, sizeof(size_t));
    size_t *c = calloc(block_size, sizeof(size_t));

    size_t *tmp_small = calloc(block_size, sizeof(size_t));
    size_t *tmp_small2 = calloc(block_size, sizeof(size_t));
    size_t *tmp_small3 = calloc(block_size, sizeof(size_t));

    size_t *intermediate_result_1 = calloc(nb_relations, sizeof(size_t));
    size_t *intermediate_result_2 = calloc(nb_relations, sizeof(size_t));
    size_t *intermediate_result_3 = calloc(nb_relations, sizeof(size_t));

    size_t *intermediate_result_4 = calloc(block_size, sizeof(size_t));

    size_t *intermediate_result_5 = calloc(nb_relations, sizeof(size_t));
    size_t *intermediate_result_6 = calloc(nb_relations, sizeof(size_t));
    size_t *intermediate_result_7 = calloc(nb_relations, sizeof(size_t));

    while (d && i <= (double)sparse_matrix.len/((double)block_size - 0.764) + 10)
    {
        multiply_sparse(sparse_matrix, nb_relations, index, V, tmp);
        sparse_multiply_transpose(sparse_matrix, tmp, Z, index);

        dense_multiply_transpose(vAv, V, Z, nb_relations, block_size);
        dense_multiply_transpose(vAAv, Z, Z, nb_relations, block_size);

        extract_columns(W_inv, &d, vAv, block_size);

        dense_multiply_transpose(tmp, V, V0, nb_relations, block_size);
        tmp2 = dense_multiply(W_inv, tmp, block_size, block_size);
        tmp = dense_multiply(V, tmp2, nb_relations, block_size);
        add_vectors(tmp2, X, tmp, nb_relations);
        memcpy(X, tmp2, nb_relations*sizeof(size_t));

        neg_d = switch_indices(d, mask);

        multiply_d(tmp_small, vAAv, d, block_size);
        multiply_d(tmp_small2, vAv, neg_d, block_size);
        add_vectors(tmp_small3, tmp_small, tmp_small2, block_size);
        c = dense_multiply(W_inv, tmp_small3, block_size, block_size);

        multiply_d(intermediate_result_1, Z, d, block_size);
        multiply_d(intermediate_result_2, V, neg_d, block_size);
        intermediate_result_3 = dense_multiply(V, c, nb_relations, block_size);
        multiply_d(intermediate_result_4, vAv, d, block_size);
        intermediate_result_5 = dense_multiply(P, intermediate_result_4, nb_relations, block_size);
        intermediate_result_6 = dense_multiply(V, W_inv, nb_relations, block_size);
        multiply_d(intermediate_result_7, P, neg_d, block_size);

        add_vectors(P, intermediate_result_6, intermediate_result_7, nb_relations);

        add_vectors(tmp, intermediate_result_1, intermediate_result_2, nb_relations);
        add_vectors(tmp2, tmp, intermediate_result_3, nb_relations);
        add_vectors(tmp, tmp2, intermediate_result_5, nb_relations);
        memcpy(V, tmp, nb_relations*sizeof(size_t));

        i++;
    }

    add_vectors(tmp, X, Y, nb_relations);
    concatenate(tmp2, tmp, V, block_size);

    mpz_t *matrix = calloc(2*block_size, sizeof(mpz_t));
    for (size_t i = 0 ; i < 2*block_size ; i++) mpz_init_set_ui(matrix[i], 0);

    mpz_t *tmp_matrix = calloc(2*block_size, sizeof(mpz_t));
    for (size_t i = 0 ; i < 2*block_size ; i++) mpz_init_set_ui(tmp_matrix[i], 0);

    multiply_sparse(sparse_matrix, nb_relations, index, tmp2, tmp);
    transpose_dense(matrix, tmp, nb_relations, 2*block_size);

    transpose_dense(tmp_matrix, tmp2, nb_relations, 2*block_size);

    solve(matrix, tmp_matrix, nb_relations, 2*block_size);

    for (size_t i = 0 ; i < 2*block_size ; i++)
    {
        if (!(mpz_cmp_ui(matrix[i], 0)))
        {
            append(output, tmp_matrix[i]);
        }
    }
    if (output->len == 0)
    {
        block_lanczos(output, sparse_matrix, nb_relations, MIN(2*block_size, 8), index);
    }
    return;
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